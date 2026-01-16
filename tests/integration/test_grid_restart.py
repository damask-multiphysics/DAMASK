# SPDX-License-Identifier: AGPL-3.0-or-later
import os
from functools import partial

import pytest
import numpy as np
import h5py

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'grid_restart'


def h5py_compare_files(iterator,asserter,
                       fname_reference,fname_current):
    with h5py.File(fname_reference,'r') as ref, \
         h5py.File(fname_current,'r') as cur:
        for (path,dset) in iterator(ref):
           if len(dset.dtype) > 0:                                                                  # cell_to
               for e in dset.dtype.fields:
                   assert np.all(dset[e]==cur[path][e])
           elif dset.dtype.char == 'S':                                                             # setup
               dset == cur[path]
           else:
               asserter(dset,cur[path])

@pytest.mark.parametrize('solver',['spectral_basic','spectral_polarization','FEM','spectral_Galerkin'])
def test_grid_restart(res_path,tmp_path,copy_files,h5py_dataset_iterator,assert_allclose,solver,petsc_version):
    grid = '27grains3x3x3' if (solver=='FEM' and petsc_version() < '3.24.1') else '8grains2x2x2'
    load = 'tensionX'
    material = 'material'

    copy_files(res_path,tmp_path)

    config_load = damask.YAML.load(res_path/f'{load}.yaml')
    config_load['solver']['mechanical'] = solver

    for mode in ['normal','restart']:
        if mode == 'restart':
            config_load['loadstep'][0]['discretization']['t']=17.
            config_load['loadstep'][0]['discretization']['N']=17
        config_load.save(tmp_path/f'{load}.yaml')

        cmd = f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {mode}'
        damask.util.run(cmd,wd=tmp_path)

        if mode == 'restart':
            config_load['loadstep'][0]['discretization']['t']=18.
            config_load['loadstep'][0]['discretization']['N']=18
            config_load.save(tmp_path/f'{load}.yaml')
            damask.util.run(cmd+' -r 17',wd=tmp_path)

    assert (damask.Table.load(tmp_path/'normal.sta').get('IterationsNeeded') ==
            damask.Table.load(tmp_path/'restart.sta').get('IterationsNeeded')).all()
    h5py_compare_files(h5py_dataset_iterator,
                       partial(assert_allclose,rtol=1.e-3),
                       tmp_path/'normal.hdf5',
                       tmp_path/'restart.hdf5')

def test_grid_restart_new_file(res_path,tmp_path,copy_files,h5py_dataset_iterator,assert_allclose):
    grid = '8grains2x2x2'
    load = 'tensionX'
    material = 'material'

    copy_files(res_path,tmp_path)

    config_load = damask.YAML.load(res_path/f'{load}.yaml')
    config_load['solver']['mechanical'] = 'spectral_Galerkin'
    config_load.save(tmp_path/f'{load}.yaml')

    damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j normal',wd=tmp_path)
    os.rename(tmp_path/'normal_restart.hdf5',tmp_path/'restart_restart.hdf5')
    damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j restart -r 17',wd=tmp_path)

    assert (damask.Table.load(tmp_path/'normal.sta').get('IterationsNeeded')[-1] ==
            damask.Table.load(tmp_path/'restart.sta').get('IterationsNeeded')[-1]).all()

    h5py_compare_files(h5py_dataset_iterator,
                       partial(assert_allclose,rtol=1.e-3),
                       tmp_path/'restart.hdf5',
                       tmp_path/'normal.hdf5')
    damask.Result(tmp_path/'restart.hdf5')


@pytest.mark.parametrize('solver',['spectral_basic','spectral_polarization','FEM','spectral_Galerkin'])
def test_grid_restart_thermo_mechanical(res_path,tmp_path,copy_files,h5py_dataset_iterator,assert_allclose,solver):
    grid = '27grains3x3x3'
    load = 'tensionX'
    material = 'material'

    copy_files(res_path,tmp_path)

    g = damask.GeomGrid.load(tmp_path/grid)
    g.initial_conditions['T'] = 300.
    g.save(tmp_path/grid)

    config_material = damask.ConfigMaterial()
    config_material['homogenization'] = {'SX': {'N_constituents': 1,'mechanical': {'type': 'pass'},
                                                                    'thermal': {'type': 'pass', 'output': ['T']}}}
    config_material = config_material.material_add(phase=['Phenopowerlaw'],
                                                   O=damask.Rotation.from_random(np.prod(g.cells),rng_seed=1),
                                                   homogenization = 'SX')
    config_material['phase']['Phenopowerlaw'] = damask.ConfigMaterial.load(res_path/'material.yaml')['phase']['Phenopowerlaw1']
    config_material['phase']['Phenopowerlaw']['rho'] = 1
    config_material['phase']['Phenopowerlaw']['thermal'] = {'C_p': 1, 'K_11': 0.0, 'K_33': 0.0,
                                                            'source': [{'type': 'externalheat', 'f': [1,1], 't': [0,100]}]}
    config_material.save(tmp_path/f'{material}.yaml')

    config_load = damask.YAML.load(res_path/f'{load}.yaml')
    config_load['solver']['mechanical'] = solver
    config_load['solver']['thermal'] = 'spectral'

    for mode in ['normal','restart']:
        if mode == 'restart':
            config_load['loadstep'][0]['discretization']['t']=17.
            config_load['loadstep'][0]['discretization']['N']=17
        config_load.save(tmp_path/f'{load}.yaml')

        cmd = f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {mode}'
        damask.util.run(cmd,wd=tmp_path)

        if mode == 'restart':
            config_load['loadstep'][0]['discretization']['t']=18.
            config_load['loadstep'][0]['discretization']['N']=18
            config_load.save(tmp_path/f'{load}.yaml')
            damask.util.run(cmd+' -r 17',wd=tmp_path)

    assert (damask.Table.load(tmp_path/'normal.sta').get('IterationsNeeded') ==
            damask.Table.load(tmp_path/'restart.sta').get('IterationsNeeded')).all()
    h5py_compare_files(h5py_dataset_iterator,
                       partial(assert_allclose,rtol=1.e-3),
                       tmp_path/'normal.hdf5',
                       tmp_path/'restart.hdf5')


def test_grid_restart_damage_simple(res_path,tmp_path,copy_files,h5py_dataset_iterator,assert_allclose):
    grid = '9x9x1_central_crack_1_layer'
    load = 'tension_brittle'
    material = 'material_brittle'

    copy_files(res_path,tmp_path, [f'{grid}.vti',f'{material}.yaml'])

    config_load = damask.YAML.load(res_path/f'{load}.yaml')

    for mode in ['normal','restart']:
        if mode == 'restart':
            config_load['loadstep'][0]['discretization']['t']=.4/9*8
            config_load['loadstep'][0]['discretization']['N']=8
        config_load.save(tmp_path/f'{load}.yaml')

        cmd = f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {mode}'
        damask.util.run(cmd,wd=tmp_path)

        if mode == 'restart':
            config_load['loadstep'][0]['discretization']['t']=.4
            config_load['loadstep'][0]['discretization']['N']=9
            config_load.save(tmp_path/f'{load}.yaml')
            damask.util.run(cmd+' -r 8',wd=tmp_path)

    assert (damask.Table.load(tmp_path/'normal.sta').get('IterationsNeeded') ==
            damask.Table.load(tmp_path/'restart.sta').get('IterationsNeeded')).all()
    h5py_compare_files(h5py_dataset_iterator,
                        partial(assert_allclose,rtol=1.e-3),
                        tmp_path/'normal.hdf5',
                        tmp_path/'restart.hdf5')


def test_grid_restart_damage_complex(res_path,tmp_path,copy_files,h5py_dataset_iterator,assert_allclose,):
    grid = 'plate_with_hole'
    load = 'isoBrittle' # 500 loadsteps
    material = 'plate_with_hole_isoBrittle_cube'

    copy_files(res_path,tmp_path, [f'{grid}.vti',f'{material}.yaml'])

    config_load = damask.YAML.load(res_path/f'{load}.yaml')

    for mode in ['normal','restart']:
        if mode == 'restart':
            config_load['loadstep'][-1]['discretization']['t'] = 0.0104
            config_load['loadstep'][-1]['discretization']['N'] = 60
            config_load['loadstep'][-1]['f_restart'] = 60 # restart record step
        config_load.save(tmp_path/f'{load}.yaml')

        cmd = f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {mode}'
        damask.util.run(cmd,wd=tmp_path)

        if mode == 'restart':
            config_load['loadstep'][-1]['discretization']['t']=0.052
            config_load['loadstep'][-1]['discretization']['N']=300
            config_load.save(tmp_path/f'{load}.yaml')
            damask.util.run(cmd+' -r 260',wd=tmp_path)

    assert (damask.Table.load(tmp_path/'normal.sta').get('IterationsNeeded') ==
            damask.Table.load(tmp_path/'restart.sta').get('IterationsNeeded')).all()
    h5py_compare_files(h5py_dataset_iterator,
                        partial(assert_allclose,rtol=1.e-3),
                        tmp_path/'normal.hdf5',
                        tmp_path/'restart.hdf5')
