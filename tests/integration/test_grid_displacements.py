# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'grid_displacements'


def test_grid_displacement(res_path,tmp_path,copy_files,assert_allclose,np_rng):

    copy_files(res_path,tmp_path)
    damask.ConfigMaterial\
        .load(tmp_path/'material.yaml') \
        .material_add(O=damask.Rotation.from_random(2,rng_seed=np_rng),phase='IF',homogenization='SX') \
        .save(tmp_path/'material.yaml')
    size = np_rng.random(3)+1.0
    cells = np_rng.integers(5,12,(3))
    damask.GeomGrid(np_rng.integers(0,2,cells),size).save(tmp_path/'random_2phase.vti')

    damask.util.run(f'damask_grid -l mixed.yaml -g random_2phase.vti -m material.yaml -w {tmp_path} '\
                     '--jobname random_2phase_mixed')

    r = damask.Result(f'{tmp_path}/random_2phase_mixed.hdf5').view(increments=-1)
    F = r.place('F').reshape(tuple(r.cells)+(3,3),order='F')
    u = r.place('u_p').reshape(tuple(r.cells)+(3,),order='F')

    assert_allclose(u,damask.grid_filters.displacement_point(r.size,F))
