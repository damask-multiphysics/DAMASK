#!/usr/bin/env python3

import os
import sys
import glob
import argparse
import shutil
from pathlib import Path
import subprocess
import shlex

sys.path.append(str(Path(__file__).resolve().parents[2]/'python/damask'))
import solver

APPLY = 'install'
RESTORE = 'uninstall'

def copy_and_patch(patch,orig,editor):
    try:
        shutil.copyfile(orig,orig.parent/patch.stem)
    except shutil.SameFileError:
        pass
    subprocess.run(shlex.split(f'patch {orig.parent/patch.stem} {patch} --backup --forward'))
    with open(orig.parent/patch.stem) as f_in:
        content = f_in.read()
    with open(orig.parent/patch.stem,'w') as f_out:
        f_out.write(content.replace('%EDITOR%',editor))


parser = argparse.ArgumentParser(
                  description=f'{APPLY.capitalize()} or {RESTORE} DAMASK modifications to MSC Marc/Mentat',
                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('command', metavar='command', nargs='?', default=APPLY, choices=[APPLY,RESTORE],
                    help=f'Mode of operation {[APPLY,RESTORE]}')
parser.add_argument('--editor', dest='editor', metavar='string', default='vi',
                    help='Name of the editor (executable) used by Marc Mentat')
parser.add_argument('--marc-root', dest='marc_root', metavar='string',
                    default=solver._marc._marc_root,
                    help='Marc root directory')
parser.add_argument('--marc-version', dest='marc_version', metavar='string',
                    default=solver._marc._marc_version,
                    help='Marc version')
parser.add_argument('--damask-root', dest='damask_root', metavar = 'string',
                    default=solver._marc._damask_root,
                    help='DAMASK root directory')

args = parser.parse_args()

damask_root  = Path(args.damask_root).expanduser()
marc_root    = Path(args.marc_root).expanduser()
marc_version = args.marc_version

matches = {'Marc_tools':  [['comp_user','comp_damask_*mp'],
                           ['run_marc','run_damask_*mp'],
                           ['include_linux64','include_linux64']],
           'Mentat_bin':  [['edit_window','edit_window'],
                           ['submit1','submit?'],
                           ['kill1','kill?']],
           'Mentat_menus':[['job_run.ms','job_run.ms']]}

for cmd in ['xvfb-run','patch'] if args.command == APPLY else ['xvfb-run'] if args.command == RESTORE else []:
    try:
        subprocess.run([cmd,'--help'], capture_output=True)
    except FileNotFoundError:
        print(f'"{cmd}" not found, please install')
        sys.exit()

if args.command == APPLY:
    print('patching files...')

    for directory in glob.glob(str(damask_root/'install/MarcMentat'/marc_version/'*')):
        for orig, mods in matches[Path(directory).name]:
            product,subfolder = (marc_root/directory).name.split('_')
            orig = marc_root/f'{product.lower()}{marc_version}/{subfolder}/{orig}'
            for patch in glob.glob(f'{directory}/{mods}.patch'):
                copy_and_patch(Path(patch),orig,args.editor)

elif args.command == RESTORE:
    print('deleting modified files...')

    for file in (glob.glob(str(marc_root/f'marc{marc_version}/tools/*_damask*')) +
                 glob.glob(str(marc_root/f'mentat{marc_version}/bin/kill[4-6]')) +
                 glob.glob(str(marc_root/f'mentat{marc_version}/bin/submit[4-6]*'))):
        os.remove(file)

    print('restoring original files...')

    for file in (glob.glob(str(marc_root/f'marc{marc_version}/tools/include_linux64.orig')) +
                 glob.glob(str(marc_root/f'mentat{marc_version}/bin/edit_window.orig')) +
                 glob.glob(str(marc_root/f'mentat{marc_version}/menus/job_run.ms.orig'))):
        shutil.copyfile(file,Path(file).with_suffix(''))
        os.remove(file)
else:
    print('skipping file modifications...')

print('compiling Mentat menu binaries...')

executable = marc_root/f'mentat{marc_version}/bin/mentat'
menu_file  = marc_root/f'mentat{marc_version}/menus/linux64/main.msb'
subprocess.run(shlex.split(f'xvfb-run -a {executable} -compile {menu_file}'))

print('setting file access rights...')

for file in (
             (glob.glob(str(marc_root/f'marc{marc_version}/tools/*_damask*')) +
              glob.glob(str(marc_root/f'mentat{marc_version}/bin/kill[4-6]')) +
              glob.glob(str(marc_root/f'mentat{marc_version}/bin/submit[4-6]')) if args.command == APPLY else []) +
             glob.glob(str(marc_root/f'mentat{marc_version}/bin/edit_window'))
            ):
    os.chmod(file , 0o755)
