# SPDX-License-Identifier: AGPL-3.0-or-later
import tempfile
import shutil
from pathlib import Path

import h5py
import numpy as np
import pytest

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'compile_mesh'

ref_file = None


@pytest.mark.parametrize('openMP',['OFF','ON'])
@pytest.mark.parametrize('optimization',['OFF','DEFENSIVE','AGGRESSIVE'])
def test_compile_mesh(damask_root,res_path,tmp_path,copy_files,h5py_dataset_iterator,assert_allclose,
                      openMP,optimization):
    global ref_file
    cmd = 'cmake '+\
            '-D CMAKE_BUILD_TYPE=DEBUG '+\
            '-D CMAKE_VERBOSE_MAKEFILE=ON '+\
            f'-D CMAKE_INSTALL_PREFIX={tmp_path} '+\
            '-D MESH=ON '+\
            f'-D OPENMP={openMP} '+\
            f'-D OPTIMIZATION={optimization} '+\
            str(damask_root)
    damask.util.run(cmd,wd=tmp_path)
    damask.util.run('make VERBOSE=1 -j4 install',wd=tmp_path)

    copy_files(res_path,tmp_path)
    load = 'tensionX'
    mesh = 'singleCrystal'
    material = 'material'
    job = f'{mesh}_{load}'
    damask.util.run(f'./bin/damask_mesh -l {load}.yaml -g {mesh}.msh -m {material}.yaml -j {job}',wd=tmp_path)

    if ref_file is None:
        print('stored reference')
        ref_file = Path(tempfile.mkdtemp())/f'{job}.hdf5'
        shutil.copyfile(tmp_path/f'{job}.hdf5',ref_file)

    cur = h5py.File(tmp_path/f'{job}.hdf5','r')
    ref = h5py.File(ref_file,'r')
    for (path,dset) in h5py_dataset_iterator(ref):
        if len(dset.dtype) > 0:                                                                     # cell_to
            for e in dset.dtype.fields:
                assert np.all(dset[e]==cur[path][e])
        elif dset.dtype.char == 'S':                                                                # setup
            dset == cur[path]
        else:
            m = max(np.amax(np.abs(dset)),np.amax(np.abs(cur[path])))
            if m < 1: continue
            assert_allclose(dset,cur[path],rtol=.0,atol=m*1e-4)
