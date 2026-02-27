# SPDX-License-Identifier: AGPL-3.0-or-later
import tempfile
import shutil
from pathlib import Path
import pytest
import h5py
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'compile_grid'


reference_results = {}

@pytest.mark.parametrize('openMP',['OFF','ON'])
@pytest.mark.parametrize('optimization',['OFF','DEFENSIVE','AGGRESSIVE'])
def test_compile_grid(damask_root,res_path,tmp_path,copy_files,h5py_dataset_iterator,assert_allclose,
                      openMP,optimization):
    global reference_results
    cmd = 'cmake '+\
            '-D CMAKE_BUILD_TYPE=DEBUG '+\
            '-D CMAKE_VERBOSE_MAKEFILE=ON '+\
            f'-D CMAKE_INSTALL_PREFIX={tmp_path} '+\
            '-D GRID=ON '+\
            f'-D OPENMP={openMP} '+\
            f'-D OPTIMIZATION={optimization} '+\
            str(damask_root)
    damask.util.run(cmd,wd=tmp_path)
    damask.util.run('make VERBOSE=1 -j4 install',wd=tmp_path)

    copy_files(res_path,tmp_path)
    load = 'tensionX'
    grid = 'test'
    material = 'material'
    job = f'{grid}_{load}'

    for solver in ['spectral_basic','spectral_polarization','FEM','spectral_Galerkin']:

        load_config = damask.YAML.load(tmp_path/f'{load}.yaml')
        load_config['solver']['mechanical'] = solver
        load_config.save(tmp_path/f'{load}.yaml')

        out, _ = damask.util.run(f'./bin/damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',
                                 wd=tmp_path)
        assert solver in out

        if solver not in reference_results:
            print('stored reference')
            reference_results[solver] = Path(tempfile.mkdtemp())/f'{job}.hdf5'
            shutil.copyfile(tmp_path/f'{job}.hdf5',reference_results[solver])

        with h5py.File(reference_results[solver],'r') as ref, \
             h5py.File(tmp_path/f'{job}.hdf5','r') as cur:
            for (path,dset) in h5py_dataset_iterator(ref):
                if len(dset.dtype) > 0:                                                             # cell_to
                    for e in dset.dtype.fields:
                        assert np.all(dset[e]==cur[path][e])
                elif dset.dtype.char == 'S':                                                        # setup
                    dset == cur[path]
                else:
                    m = max(np.amax(np.abs(dset)),np.amax(np.abs(cur[path])))
                    if m < 1: continue
                    assert_allclose(dset,cur[path],rtol=.0,atol=m*1e-4)
