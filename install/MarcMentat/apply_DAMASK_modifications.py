#!/usr/bin/env python3

import os
import glob
import argparse
import shutil
from pathlib import Path

import damask

def copy_and_patch(patch,orig,editor):
    try:
        shutil.copyfile(orig,orig.parent/patch.stem)
    except shutil.SameFileError:
        pass
    damask.util.run(f'patch {orig.parent/patch.stem} {patch} -b')
    with open(orig.parent/patch.stem) as f_in:
        content = f_in.read()
    with open(orig.parent/patch.stem,'w') as f_out:
        f_out.write(content.replace('%EDITOR%',editor))


parser = argparse.ArgumentParser(
                  description='Apply DAMASK modification to MSC Marc/Mentat',
                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--editor', dest='editor', metavar='string', default='vi',
                    help='Name of the editor for Marc Mentat (executable)')
parser.add_argument('--marc-root', dest='marc_root', metavar='string',
                    default=damask.solver._marc._marc_root,
                    help='Marc root directory')
parser.add_argument('--marc-version', dest='marc_version', metavar='string',
                    default=damask.solver._marc._marc_version,
                    help='Marc version')
parser.add_argument('--damask-root', dest='damask_root', metavar = 'string',
                    default=damask.solver._marc._damask_root,
                    help='DAMASK root directory')

args = parser.parse_args()
marc_root = Path(args.marc_root).expanduser()
damask_root = Path(args.damask_root).expanduser()
marc_version = args.marc_version

matches = {'Marc_tools':  [['comp_user','comp_damask_*mp'],
                           ['run_marc','run_damask_*mp'],
                           ['include_linux64','include_linux64']],
           'Mentat_bin':  [['edit_window','edit_window'],
                           ['submit1','submit?'],
                           ['kill1','kill?']],
           'Mentat_menus':[['job_run.ms','job_run.ms']]}


print('patching files...')

for directory in glob.glob(str(damask_root/'install/MarcMentat'/marc_version/'*')):
    for orig, mods in matches[Path(directory).name]:
        product,subfolder = (marc_root/Path(directory)).name.split('_')
        orig = marc_root/f'{product.lower()}{marc_version}/{subfolder}/{orig}'
        for patch in glob.glob(f'{directory}/{mods}.patch'):
            copy_and_patch(Path(patch),orig,args.editor)

print('compiling Mentat menu binaries...')

executable = marc_root/f'mentat{marc_version}/bin/mentat'
menu_file  = marc_root/f'mentat{marc_version}/menus/linux64/main.msb'
os.system(f'xvfb-run -a {executable} -compile {menu_file}')

print('setting file access rights...')

for file in (glob.glob(str(marc_root/f'marc{marc_version}/tools/*_damask*')) +
             glob.glob(str(marc_root/f'mentat{marc_version}/bin/kill[4-6]')) +
             glob.glob(str(marc_root/f'mentat{marc_version}/bin/submit[4-6]'))):
    os.chmod(file , 0o755)
