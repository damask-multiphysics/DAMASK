# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'damage_scaling'


def test_resolution(res_path,tmp_path,copy_files):
    copy_files(res_path,tmp_path)
    for grid in ['40','80']:
        load = f'tension{grid}'
        material = 'material'
        job = f'{grid}_{load}'
        damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job} -n numerics.yaml',
                        wd=tmp_path)
        r = damask.Result(f'{tmp_path}/{job}.hdf5')
        P = [np.average(P_[:,0,0]) for P_ in r.place('P').values()]
        assert np.max(P) > 20e6
        assert P[-1] < 17e6
