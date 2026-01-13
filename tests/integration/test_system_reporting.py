# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'system_reporting'


def test_system_reporting(copy_files,res_path,tmp_path,np_rng):
    grid = 'simple'
    load = 'shearZX'
    material = 'material'

    copy_files(res_path,tmp_path)

    N_sl = np.rint(np_rng.random(5)*np.array([3,3,6,12,6])).astype(int)
    N_tw = np.rint(np_rng.random(4)*np.array([6,6,6,6])).astype(int)

    config = damask.ConfigMaterial.load(tmp_path/f'{material}.yaml')
    config['phase']['Mg']['mechanical']['plastic']['N_sl'] = N_sl.tolist()
    config['phase']['Mg']['mechanical']['plastic']['N_tw'] = N_tw.tolist()
    config['material'][0]['constituents'][0]['O']=damask.Rotation.from_random(rng_seed=np_rng)
    config.save(tmp_path/f'{material}.yaml')

    damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {grid}_{load}',
                    wd=tmp_path)

    r = damask.Result(tmp_path/f'{grid}_{load}.hdf5')
    assert N_sl.sum() == 0 or len(r.view(increments=-1).get('xi_sl').dtype.metadata['systems']) == N_sl.sum()
    assert N_tw.sum() == 0 or len(r.view(increments=-1).get('xi_tw').dtype.metadata['systems']) == N_tw.sum()
