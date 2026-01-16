# SPDX-License-Identifier: AGPL-3.0-or-later
import os

import pytest
import h5py
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'Marc_element_lib'


@pytest.mark.parametrize('elem',[6,7,11,21,27,54,57,117,125,127,134,136])
def test_marc_element_lib(damask_root,res_path,tmp_path,copy_files,assert_allclose,update,elem):
    job = 'planeStrain' if elem in [125,6,27,54,57,11] else 'tension'
    copy_files(res_path,tmp_path,[f'elemtype{elem}_{job}.dat','material.yaml'])
    s = damask.solver.Marc(damask_root=damask_root)
    os.chdir(tmp_path)
    s.submit_job(model=f'elemtype{elem}', job=job)

    if update:
        with h5py.File(res_path/'Reference.hdf5','r+') as ref, \
             h5py.File(tmp_path/f'elemtype{elem}_{job}.hdf5','r') as cur:
            for dataset in ['T_e', 'T_c','x_n','x_p']:
                loc_cur = f'/geometry/{dataset}'
                loc_ref = f'/geometry_{elem}/{dataset}'
                ref[loc_ref][()]=cur[loc_cur][()]
            for displacement in ['u_n','u_p']:
                last_inc=[s for s in list(ref.keys()) if f'lastInc_{elem}' in s][0]
                loc_cur = f'increment_{last_inc[3:6]}/geometry/{displacement}'
                loc_ref= last_inc+f'/geometry/{displacement}'
                ref[loc_ref][()]=cur[loc_cur][()]

    with h5py.File(res_path/'Reference.hdf5','r') as ref, \
         h5py.File(tmp_path/f'elemtype{elem}_{job}.hdf5','r') as cur:
        for dataset in ['T_e', 'T_c','x_n','x_p']:
            loc_cur = f'/geometry/{dataset}'
            loc_ref = f'/geometry_{elem}/{dataset}'
            assert np.all(ref[loc_ref][()] == cur[loc_cur][()])
        for displacement in ['u_n','u_p']:
            last_inc=[s for s in list(ref.keys()) if f'lastInc_{elem}' in s][0]
            loc_cur = f'increment_{last_inc[3:6]}/geometry/{displacement}'
            loc_ref= last_inc+f'/geometry/{displacement}'
            assert_allclose(ref[loc_ref][()],cur[loc_cur][()])
