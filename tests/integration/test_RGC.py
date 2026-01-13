# SPDX-License-Identifier: AGPL-3.0-or-later
import shutil

import pytest
import numpy as np
import h5py

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'RGC'


def test_RGC_detect_changes(res_path,tmp_path,copy_files,h5py_dataset_iterator,assert_allclose,update):
    copy_files(res_path,tmp_path)
    load = 'tensionX'
    grid = '3grains3x3x3'
    material = 'material'
    job = f'{grid}_{load}'

    stdout,stderr = damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',
                                    wd=tmp_path)
    print(f'{stdout}\n{stderr}')

    if update:
        shutil.copy(tmp_path/f'{job}.hdf5',res_path)
    with h5py.File(res_path/f'{job}.hdf5','r') as ref, \
         h5py.File(tmp_path/f'{job}.hdf5','r') as cur:
        for (path,dset) in h5py_dataset_iterator(ref):
            if path.startswith('setup'): continue
            if len(dset.dtype) > 0:
                print(f'comparing {path}/{dset} (compound data type)')
                for e in dset.dtype.fields:
                    assert np.all(dset[e]==cur[path][e])
            else:
                print(f'comparing {path}/{dset}')
                atol=max(1e-8*max(np.amax(np.abs(dset)),np.amax(np.abs(cur[path]))),1e-12)
                assert_allclose(dset,cur[path],rtol=.0,atol=atol)
