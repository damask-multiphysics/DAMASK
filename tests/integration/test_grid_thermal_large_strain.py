# SPDX-License-Identifier: AGPL-3.0-or-later
import os

import pytest
import numpy as np
import yaml
import h5py

import damask

"""
Test correctness of large strain formulation.

A (10,10,10) grid of initial size (1,1,1) and the same grid
with initial size (2,1,.5) deformed to size (1,1,1) should
give the same behavior.

"""

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'grid_thermal_large_strain'


def test_grid_thermal_large_strain(res_path,tmp_path,copy_files,assert_allclose,np_rng):

    grid = [10]*3

    os.chdir(tmp_path)
    copy_files(res_path,tmp_path)

    g = damask.GeomGrid(np.zeros(grid),np.ones(3)*1e-4)
    T = np_rng.uniform(100.0,500.0,size=grid)
    g.initial_conditions['T'] = T
    g.save('equispaced')

    load = damask.YAML(
          '''
          solver:
            mechanical: spectral_basic
            thermal: spectral
          loadstep:
            - boundary_conditions:
                mechanical:
                  dot_F: [[0, 0, 0],
                          [0, 0, 0],
                          [0, 0, 0]]
              discretization:
                t: 1
                N: 50
                r: 1.02
            ''')
    load.save('none.yaml')


    damask.util.run('DAMASK_grid -l none.yaml -g equispaced.vti -m material.yaml -j equispaced_none')

    compress = yaml.safe_load(
              '''
              boundary_conditions:
                mechanical:
                  F: [[.5, 0, 0],
                      [0, 1, 0],
                      [0, 0, 2]]
              discretization:
                t: 1
                N: 10
              f_restart: 10
              f_out: none
              ''')


    load['loadstep'].insert(0,compress)
    load.save('pre-compress.yaml')
    g.size *= (2,1,.5)
    g.save('distorted')
    damask.util.run('DAMASK_grid -l pre-compress.yaml -g distorted.vti -m material.yaml '\
                    '-j distorted_pre-compress')

    with h5py.File('distorted_pre-compress_restart.hdf5','a') as f:
        f['solver']['T'][...] = T.reshape(-1,1,order='F')
        f['solver']['T_lastinc'][...] = T.reshape(-1,1,order='F')

    with h5py.File('distorted_pre-compress.hdf5','a') as f:
        for inc in range(11,61):
            del(f[f'increment_{inc}'])

    damask.util.run('DAMASK_grid -l pre-compress.yaml -g distorted.vti -m material.yaml '\
                    '-j distorted_pre-compress --restart 10')

    r_equispaced = damask.Result('equispaced_none.hdf5')
    r_distorted = damask.Result('distorted_pre-compress.hdf5')
    for i in range(1,51):
        assert_allclose(r_equispaced.view(increments=i).get('T'),
                        r_distorted.view(increments=i+10).get('T'))
