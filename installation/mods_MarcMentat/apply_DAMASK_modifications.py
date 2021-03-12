#!/usr/bin/env python3

import os
import glob
import argparse
from pathlib import Path

msc_version = os.environ['MSC_VERSION']
msc_root    = Path(os.environ['MSC_ROOT'])
damask_root = Path(os.environ['DAMASK_ROOT'])

parser = argparse.ArgumentParser(
     description='Apply DAMASK modification to MSC.Marc/Mentat',
     epilog = f'MSC_ROOT={msc_root} and MSC_VERSION={msc_version} (from {damask_root}/env/CONFIG)')
parser.add_argument('--editor', dest='editor', metavar='string', default='vi',
                    help='Name of the editor for MSC.Mentat (executable)')


def copy_and_replace(in_file,dst):
    with open(in_file) as f:
        content = f.read()
    content = content.replace('%INSTALLDIR%',str(msc_root))
    content = content.replace('%VERSION%', msc_version)
    content = content.replace('%EDITOR%', parser.parse_args().editor)
    with open(dst/Path(in_file).name,'w') as f:
        f.write(content)


print('adapting Marc tools...\n')

src = damask_root/f'installation/mods_MarcMentat/{msc_version}/Marc_tools'
dst = msc_root/f'marc{msc_version}/tools'
for in_file in glob.glob(str(src/'*damask*')) + [str(src/'include_linux64')]:
    copy_and_replace(in_file,dst)


print('adapting Mentat scripts and menus...\n')

src = damask_root/f'installation/mods_MarcMentat/{msc_version}/Mentat_bin'
dst = msc_root/f'mentat{msc_version}/bin'
for in_file in glob.glob(str(src/'*[!.original]')):
    copy_and_replace(in_file,dst)

src = damask_root/f'installation/mods_MarcMentat/{msc_version}/Mentat_menus'
dst = msc_root/f'mentat{msc_version}/menus'
for in_file in glob.glob(str(src/'job_run.ms')):
    copy_and_replace(in_file,dst)


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
