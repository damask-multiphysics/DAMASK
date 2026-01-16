# SPDX-License-Identifier: AGPL-3.0-or-later
import os

import pytest
import h5py
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'Marc_ddm'


def test_marc_ddm(damask_root,res_path,tmp_path,copy_files,assert_allclose,update):
    copy_files(res_path,tmp_path,['r-wert_new_R-Wert.dat','1r-wert_new_R-Wert.dat',\
                                  '2r-wert_new_R-Wert.dat','3r-wert_new_R-Wert.dat','material.yaml'])
    s = damask.solver.Marc(damask_root=damask_root)
    os.chdir(tmp_path)
    s.submit_job(model='r-wert_new', job='R-Wert', domains=3)

    for domain in ['1','2','3']:
        if update:
            with h5py.File(res_path/'Reference.hdf5','r+') as ref, \
                 h5py.File(tmp_path/f'{domain}r-wert_new_R-Wert.hdf5','r') as cur:
                for dataset in ['T_e', 'T_c','x_n','x_p']:
                    loc_cur = f'/geometry/{dataset}'
                    loc_ref = f'/geometry_{domain}/{dataset}'
                    ref[loc_ref][()]=cur[loc_cur][()]
                for displacement in ['u_n','u_p']:
                    last_inc=[s for s in list(ref.keys()) if f'lastInc_{domain}' in s][0]
                    loc_cur = f'increment_{last_inc[3:5]}/geometry/{displacement}'
                    loc_ref= last_inc+f'/geometry/{displacement}'
                    ref[loc_ref][()]=cur[loc_cur][()]

        with h5py.File(res_path/'Reference.hdf5','r') as ref, \
             h5py.File(tmp_path/f'{domain}r-wert_new_R-Wert.hdf5','r') as cur:
            for dataset in ['T_e', 'T_c','x_n','x_p']:
                loc_cur = f'/geometry/{dataset}'
                loc_ref = f'/geometry_{domain}/{dataset}'
                assert np.all(ref[loc_ref][()] == cur[loc_cur][()])
            for displacement in ['u_n','u_p']:
                last_inc=[s for s in list(ref.keys()) if f'lastInc_{domain}' in s][0]
                loc_cur = f'increment_{last_inc[3:5]}/geometry/{displacement}'
                loc_ref= last_inc+f'/geometry/{displacement}'
                assert_allclose(ref[loc_ref][()],cur[loc_cur][()])
