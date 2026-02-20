# SPDX-License-Identifier: AGPL-3.0-or-later
import os
import shutil
from functools import partial

import pytest
import numpy as np
import h5py

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'damage'


def h5py_compare_files(iterator,asserter,
                       fname_reference,fname_current):
    ref=h5py.File(fname_reference,'r')
    cur=h5py.File(fname_current,'r')
    for (path,dset) in iterator(ref):
       if path.startswith('setup'): continue
       if len(dset.dtype) > 0:
           for e in dset.dtype.fields:
               assert np.all(dset[e]==cur[path][e])
       else:
           if os.path.split(path)[1] == 'P':
               continue
           asserter(dset,cur[path])

@pytest.mark.parametrize('ori',['cube','rotCube'])
@pytest.mark.parametrize('struct',['fcc','bcc'])
def test_aniso_brittle(res_path,tmp_path,copy_files,h5py_dataset_iterator,assert_allclose,update,struct,ori):
    grid = 'plate_with_hole'
    load = f'anisoBrittle_{ori}'
    material = 'material'
    job = f'{grid}_{load}'

    copy_files(res_path,tmp_path,[f'{load}.yaml',f'{grid}.vti'])
    mat_config = damask.ConfigMaterial.load(res_path/f'{grid}_anisoBrittle_{struct}_{ori}.yaml')
    mat_config.save(tmp_path/f'{material}.yaml')

    damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',
                    wd=tmp_path)

    if update:
        shutil.copyfile(tmp_path/f'{job}.hdf5',res_path/f'{job}.hdf5')

    h5py_compare_files(h5py_dataset_iterator,
                       partial(assert_allclose,atol=0.0075),
                       res_path/f'{job}.hdf5',
                       tmp_path/f'{job}.hdf5')


@pytest.mark.parametrize('ori',[('cube'),('rotCube')])
def test_iso_brittle(res_path,tmp_path,copy_files,h5py_dataset_iterator,assert_allclose,update,ori):
    grid = 'plate_with_hole'
    load = 'isoBrittle'
    material = 'material'
    job = f'{grid}_{load}'

    copy_files(res_path,tmp_path,[f'{load}.yaml',f'{grid}.vti'])
    mat_config = damask.ConfigMaterial.load(res_path/f'{job}_{ori}.yaml')
    mat_config.save(tmp_path/f'{material}.yaml')

    damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',
                    wd=tmp_path)
    if update:
        shutil.copyfile(tmp_path/f'{job}.hdf5',res_path/f'{job}.hdf5')

    h5py_compare_files(h5py_dataset_iterator,
                       partial(assert_allclose,atol=0.005),
                       res_path/f'{job}.hdf5',
                       tmp_path/f'{job}.hdf5')
