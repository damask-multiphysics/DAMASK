# SPDX-License-Identifier: AGPL-3.0-or-later
import os
import shutil

import pytest

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'grid_CLI_flags'


@pytest.mark.parametrize('cmd',\
    [' --loadcase tensionX.load             --geometry test.vti              --materialconfig mat.yaml',
     ' -l=load/tensionX.yaml                 -g=geometry/tst.vti             -m=material/matconf.yaml',
     ' -l nowd/tensionX.yaml                 -g nowd/test.vti                -m nowd/mat.yaml',
     ' --load <path>/t.load                 --geom <path>/t.vti              --material <path>/m.yaml',
     ' --loadcase nowd1/nowd2/tensionX.load --geometry nowd1/nowd2/test.vti  --material=nowd1/nowd2/mat.yaml',
     ' --loadcase nowd/tens.load            --geometry=nowd3/test.vti        --material nowd3/mat.yaml --workingdirectory wd1/wd2/',
     ' -l=nwd/tension.l                      -g nwd/test.vti                 -m nwd/mat.yaml           -w <path>/x',
     ' -l ../tnsX.l                          -g nwd/test.v                   -m nwd/mat.yaml           --workingdir <path>/wd1',
     ' -l directory/tensionX.yaml            -g ../t.v                       -m=../mat.m               -w <path>/wd1/wd2',
     ' -l=../tensionX.l                      -g=../../test.v                 -m=../../mat.m            -w=wd1/wd2'])

def test_grid_CLI_flags(res_path,tmp_path,cmd):
    cmd = cmd.replace('<path>',str(tmp_path))

    alt_keys = { 'loadcase'   :[' -l',' --load',' --loadcase'],
                 'geometry'   :[' -g',' --geom',' --geometry'],
                 'material'   :[' -m',' --material',' --materialconfig'],
                 'workingdir' :[' -w',' --wd',' --workingdir',' --workingdirectory'] }
    args = {}
    for k,v in alt_keys.items():
        for a in v:
            if a in cmd: args[k] = cmd.replace('=',' ').split(a)[-1].strip().split()[0]
    if 'workingdir' not in args:
        cwd = str(tmp_path)
    else:
        cwd = str(tmp_path/args['workingdir'])
        os.makedirs(cwd)

    load_path, load_file = os.path.split(os.path.normpath(os.path.join(cwd,args['loadcase'])))
    os.makedirs(load_path,exist_ok=True)
    geom_path, geom_file = os.path.split(os.path.normpath(os.path.join(cwd,args['geometry'])))
    os.makedirs(geom_path,exist_ok=True)
    material_path, material_file = os.path.split(os.path.normpath(os.path.join(cwd,args['material'])))
    os.makedirs(material_path,exist_ok=True)

    job = os.path.splitext(geom_file)[0]+'_'+os.path.splitext(load_file)[0]

    shutil.copy2(res_path/'4x4x4homog.vti',geom_path+'/'+geom_file)
    shutil.copy2(res_path/'tensionX.yaml', load_path+'/'+load_file)
    shutil.copy2(res_path/'material.yaml', material_path+'/'+material_file)
    damask.util.run(f'DAMASK_grid {cmd} -j {job}',wd=tmp_path)

    assert os.path.isfile(os.path.join(cwd,f'{job}.hdf5'))
