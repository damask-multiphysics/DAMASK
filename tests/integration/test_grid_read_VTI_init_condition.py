# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'grid_read_VTI_init_condition'


def test_grid_read_VTI_init_phi(res_path,tmp_path,copy_files,assert_allclose,np_rng):
    load = 'tension'
    grid = 'cube3x3x3'
    material = 'material'
    job = f'{grid}_{load}'

    g = damask.GeomGrid(np.zeros([3,3,3]), [3,3,3])
    phi_0 = np.ones_like(g.material,dtype=float)
    phi_0.ravel()[phi_0.size//2] = np_rng.random()                                                  # random damage at central voxel
    g.initial_conditions['phi'] = phi_0
    g.save(tmp_path/f'{grid}')

    copy_files(res_path,tmp_path)

    damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',
                    wd=tmp_path)

    result = damask.Result(tmp_path/f'{job}.hdf5')
    result_phi_0 = result.place('phi')['increment_0'].data

    assert_allclose(phi_0.ravel(), result_phi_0)
