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
    return res_path_base/'grid_mech_consistency'


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
           if os.path.split(path)[1] in  ['u_n', 'u_p']:
               continue
           asserter(dset,cur[path])

# this class attribute is not reset when using pytest-repeat
ref = {}

@pytest.mark.parametrize('solver',['spectral_basic','spectral_polarization','spectral_Galerkin','FEM'])
@pytest.mark.parametrize('grid, load, material', [('g_333', 'tensionY', 'material_elast'),
                                                  ('g_444', 'tensionX', 'material_plast'),])
def test_grid_mech_consistency(res_path,tmp_path,copy_files,h5py_dataset_iterator,assert_allclose,
                               solver,grid,load,material):
    job = f'{grid}_{load}_{material}'

    copy_files(res_path,tmp_path,[f'{grid}.vti',f'{material}.yaml'])
    l = damask.LoadcaseGrid.load(res_path/f'{load}.yaml')
    l['solver']['mechanical'] = solver
    l.save(tmp_path/f'{load}.yaml')

    damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml',
                    wd=tmp_path)

    if job in ref:
        h5py_compare_files(h5py_dataset_iterator,
                           partial(assert_allclose,atol=1e-2),
                           ref[job],
                           tmp_path/f'{job}.hdf5')
    else:
        ref[job] = tmp_path/f'{job}.hdf5'
