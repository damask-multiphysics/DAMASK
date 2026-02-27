# SPDX-License-Identifier: AGPL-3.0-or-later
import os
env = os.environ                                                                                    # without MPI (might come with h5py)

import pytest
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'mesh_MPI'


@pytest.mark.xfail # some combinations for p_s==2 failing with petsc\
@pytest.mark.parametrize('n',[2,3])
def test_mesh_MPI(res_path,tmp_path,copy_files,mpi_launcher,np_rng,n):

    load = 'tensionX'
    mesh = 'bicrystal_3D'
    material = 'material'
    job = f'{mesh}_{load}'

    p_s = np_rng.integers(1,6)
    p_i = np_rng.integers(max(1,p_s-1),p_s+1)

    for mode in ['serial','parallel']:
        cwd = tmp_path/mode

        os.mkdir(cwd)
        copy_files(res_path,cwd)
        numerics_config = damask.YAML.load(res_path/'numerics.yaml')
        numerics_config['solver']['mesh']['p_i'] = p_i
        numerics_config['solver']['mesh']['p_s'] = p_s
        numerics_config.save(cwd/'numerics.yaml')

        damask.util.run((f'{mpi_launcher} -n {n} ' if mode == 'parallel' else '') +
                         f'damask_mesh -l {load}.yaml -g {mesh}.msh -m {material}.yaml '+
                         f'-n numerics.yaml -j {job}',wd=cwd,env=env)
        damask.Result(tmp_path/f'{mode}/{job}.hdf5').export_VTK('*',mode='point')
