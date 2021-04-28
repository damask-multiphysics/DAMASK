#!/usr/bin/env python3

import os
import glob
import argparse
from pathlib import Path

import damask

def copy_and_replace(in_file,dst,msc_root,editor):
    with open(in_file) as f_in, open(dst/Path(in_file).name,'w') as f_out:
        f_out.write(f_in.read().replace('%INSTALLDIR%',msc_root).replace('%EDITOR%',editor))


parser = argparse.ArgumentParser(
                  description='Apply DAMASK modification to MSC.Marc/Mentat',
                  prog = Path(__file__).name,
                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--editor', dest='editor', metavar='string', default='vi',
                    help='Name of the editor for MSC.Mentat (executable)')
parser.add_argument('--msc-root', dest='msc_root', metavar='string',
                    default=damask.solver._marc._msc_root,
                    help='MSC.Marc/Mentat root directory')
parser.add_argument('--msc-version', dest='msc_version', type=float, metavar='float',
                    default=damask.solver._marc._msc_version,
                    help='MSC.Marc/Mentat version')
parser.add_argument('--damask-root', dest='damask_root', metavar = 'string',
                    default=damask.solver._marc._damask_root,
                    help='DAMASK root directory')

args = parser.parse_args()
msc_root    = Path(args.msc_root)
damask_root = Path(args.damask_root)
msc_version = args.msc_version


print('adapting Marc tools...\n')

src = damask_root/f'installation/mods_MarcMentat/{msc_version}/Marc_tools'
dst = msc_root/f'marc{msc_version}/tools'
for in_file in glob.glob(str(src/'*damask*')) + [str(src/'include_linux64')]:
    copy_and_replace(in_file,dst,args.msc_root,args.editor)


print('adapting Mentat scripts and menus...\n')

src = damask_root/f'installation/mods_MarcMentat/{msc_version}/Mentat_bin'
dst = msc_root/f'mentat{msc_version}/bin'
for in_file in glob.glob(str(src/'*[!.original]')):
    copy_and_replace(in_file,dst,args.msc_root,args.editor)

src = damask_root/f'installation/mods_MarcMentat/{msc_version}/Mentat_menus'
dst = msc_root/f'mentat{msc_version}/menus'
for in_file in glob.glob(str(src/'job_run.ms')):
    copy_and_replace(in_file,dst,args.msc_root,args.editor)


print('compiling Mentat menu binaries...')

executable = msc_root/f'mentat{msc_version}/bin/mentat'
menu_file  = msc_root/f'mentat{msc_version}/menus/linux64/main.msb'
os.system(f'xvfb-run {executable} -compile {menu_file}')


print('setting file access rights...\n')

for pattern in [msc_root/f'marc{msc_version}/tools/*damask*',
                msc_root/f'mentat{msc_version}/bin/submit?',
                msc_root/f'mentat{msc_version}/bin/kill?']:
    for f in glob.glob(str(pattern)):
        os.chmod(f,0o755)
