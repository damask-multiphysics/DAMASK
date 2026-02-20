# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'thermal_isotropic'


def test_thermal_isotropic(res_path,tmp_path,copy_files,np_rng):
    grid = 'inclusion'
    load = 'no_deformation'
    material = 'material'
    copy_files(res_path,tmp_path)

    g = damask.GeomGrid.load(tmp_path/f'{grid}.vti')
    g.initial_conditions['T']=np.where(g.material==0,np_rng.random()*400,np_rng.random()*400)+200.
    g.save(tmp_path/f'{grid}.vti')

    damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {grid}_{load} -n numerics.yaml',
                    wd=tmp_path)
    r_cube = damask.Result(tmp_path/f'{grid}_{load}.hdf5').place(['T','O']).values()

    m = damask.ConfigMaterial.load(tmp_path/f'{material}.yaml')
    m['material'][0]['constituents'][0]['O'] = damask.Rotation.from_random(rng_seed=np_rng)
    m['material'][1]['constituents'][0]['O'] = damask.Rotation.from_random(rng_seed=np_rng)
    m.save(tmp_path/f'{material}.yaml')

    damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {grid}_{load} -n numerics.yaml',
                    wd=tmp_path)
    r_random = damask.Result(tmp_path/f'{grid}_{load}.hdf5').place(['T','O']).values()

    for cube,random in zip(r_cube,r_random):
        assert (        np.allclose(cube['thermal'],random['thermal'])
                and not np.allclose(cube['mechanical'],random['mechanical']))
