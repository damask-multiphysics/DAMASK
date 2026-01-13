# SPDX-License-Identifier: AGPL-3.0-or-later
import os
import shutil

import pytest
import h5py
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'plastic_detect_changes'


def h5py_dataset_iterator(g, prefix=''):
    for key,item in g.items():
        path = os.path.join(prefix, key)
        if isinstance(item, h5py.Dataset): # test for dataset
            yield (path, item)
        elif isinstance(item, h5py.Group): # test for group (go down)
            yield from h5py_dataset_iterator(item, path)


def h5py_compare_files(fname_reference,fname_current,assert_allclose):
    with h5py.File(fname_reference,'r') as ref, \
         h5py.File(fname_current,'r') as cur:
        for (path,dset) in h5py_dataset_iterator(ref):
            if path.startswith('setup'): continue
            if len(dset.dtype) > 0:
                for e in dset.dtype.fields:
                    assert np.all(dset[e]==cur[path][e])
            elif dset.dtype == 'int32':
                assert dset.size - np.count_nonzero(cur[path][()]==ref[path][()]) <= dset.size*0.0001
            else:
                if ref[path].size < 1: continue
                rtol = 0.025
                atol = max(np.abs(dset).max()*1e-5,1e-6)
                assert_allclose(cur[path],ref[path],rtol=rtol,atol=atol)

@pytest.mark.parametrize('model',[
                                  'dislotungsten','dislotungsten_nonSchmid',
                                  'dislotwin',
                                  'kinehardening',
                                  'nonlocalbcc', 'nonlocalfcc',
                                  'phenopowerlawbcc', 'phenopowerlawhex',
                                 ])
def test_plastic_detect_changes(tmp_path,res_path,copy_files,assert_allclose,update,model):

    load = 'tensionX'
    grid = '17grains'
    material = 'material'
    job = f'{grid}_{load}'

    copy_files(res_path,tmp_path)
    m = damask.ConfigMaterial.load(tmp_path/f'{material}.yaml')
    m.material_rename_phase({'n/a':model}).save(tmp_path/f'{material}.yaml')

    stdout,stderr = damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -n numerics.yaml -j {job}',
                                    wd=tmp_path)

    if update:
        shutil.copyfile(tmp_path/f'{job}.hdf5',
                        res_path/f'{job}_{model}.hdf5')

    h5py_compare_files(res_path/f'{job}_{model}.hdf5',
                        tmp_path/f'{job}.hdf5',
                        assert_allclose)
