# SPDX-License-Identifier: AGPL-3.0-or-later
import os

import pytest
import h5py

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'Marc_output'

def test_marc_output(damask_root,res_path,tmp_path,copy_files,assert_allclose,update):
    copy_files(res_path,tmp_path,['r-wert_new_R-Wert.dat','1r-wert_new_R-Wert.dat',\
                                  '2r-wert_new_R-Wert.dat','3r-wert_new_R-Wert.dat','material.yaml'])
    s = damask.solver.Marc(damask_root=damask_root)
    os.chdir(tmp_path)
    s.submit_job(model='r-wert_new', job='R-Wert', domains=3)

    for domain in ['1','2','3']:
        f = h5py.File(tmp_path/f'{domain}r-wert_new_R-Wert.hdf5','r')
        for inc in ['2','4','6','8','10','13','16','19','20']:
            group = f'increment_{inc}'
            assert f[group]
        for inc in ['1','3','5','7','9','11','12','14','15','17','18']:
            group = f'increment_{inc}'
            assert not group in f
