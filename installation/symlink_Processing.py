#!/usr/bin/env python3

# Makes postprocessing routines accessible from everywhere.
import sys
from pathlib import Path
import os

bin_dir = Path(os.environ['DAMASK_ROOT'])/'bin'

if not bin_dir.exists():
    bin_dir.mkdir()


sys.stdout.write('\nsymbolic linking...\n')
for sub_dir in ['pre','post']:
    the_dir = Path(os.environ['DAMASK_ROOT'])/'processing'/sub_dir

    for the_file in the_dir.glob('*.py'):
        src = the_dir/the_file
        dst = bin_dir/Path(the_file.with_suffix('').name)
        if dst.is_file(): dst.unlink() # dst.unlink(True) for Python >3.8
        dst.symlink_to(src)


sys.stdout.write('\npruning broken links...\n')
for filename in bin_dir.glob('*'):
    if not filename.is_file():
        filename.unlink()
