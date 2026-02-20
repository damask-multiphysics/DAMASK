# SPDX-License-Identifier: AGPL-3.0-or-later
import os
env = os.environ                                                                                    # without MPI (might come with h5py)

import pytest
import numpy as np
import h5py

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'grid_restart_MPI'


def h5py_compare_files_statistically(iterator,asserter,
                                     fname_reference,fname_current,derivative):
    with h5py.File(fname_reference,'r') as ref, \
         h5py.File(fname_current,'r') as cur:
        for (path,dset) in iterator(ref):
           if len(dset.dtype) > 0:                                                                  # cell_to
               for e in dset.dtype.fields:
                   assert np.all(dset[e]==cur[path][e])
           elif dset.dtype.char == 'S':                                                             # setup
               dset == cur[path]
           else:
               if derivative is not None and 'geometry' in path: continue                           # coordinates are wrong/erratic
               if np.prod(dset.shape) == 0: continue

               asserter(np.mean(dset),np.mean(cur[path]))
               asserter(np.std(dset),np.std(cur[path]))

@pytest.mark.parametrize('solver,derivative',[('spectral_basic','continuous'),
                                              ('spectral_basic','FWBW_difference'),
                                              ('spectral_basic','central_difference'),
                                              ('spectral_polarization','continuous'),
                                              ('spectral_polarization','FWBW_difference'),
                                              ('spectral_polarization','central_difference'),
                                              ('spectral_Galerkin','continuous'),
                                              ('FEM',None)])
@pytest.mark.parametrize('N_proc',[1,2])
def test_grid_restart_MPI(res_path,tmp_path,copy_files,mpi_launcher,h5py_dataset_iterator,assert_allclose,
                          solver,derivative,N_proc):
    grid = 'test.vti'
    load = 'tensionX.yaml'
    material = 'material.yaml'
    numerics = 'numerics.yaml'

    copy_files(res_path,tmp_path)

    num_config = damask.YAML.load(res_path/'numerics.yaml')
    if derivative is not None: num_config['solver']['grid']['FFT']['derivative'] = derivative
    num_config.save(tmp_path/'numerics.yaml')

    load_config = damask.YAML.load(res_path/load)
    load_config['solver']['mechanical'] = solver
    load_config.save(tmp_path/load)

    for mode in ['normal','restart']:

        cmd = f'{mpi_launcher} -n {N_proc} damask_grid -l {load} -g {grid} -m {material} '\
                f'-n {numerics} -j {mode} --wd {str(tmp_path)}'
        out, _ = damask.util.run(cmd,env=env)
        assert solver in out

        if mode == 'restart':
            with h5py.File(tmp_path/f'{mode}.hdf5','r+') as f:
                del f['increment_25']
            t = damask.Table.load(tmp_path/f'{mode}.sta')
            damask.Table(t.shapes,t.data.to_numpy()[:17]).save(tmp_path/f'{mode}.sta')
            damask.util.run(cmd+' -r 17',env=env)

    assert (damask.Table.load(tmp_path/'normal.sta').get('IterationsNeeded') ==
            damask.Table.load(tmp_path/'restart.sta').get('IterationsNeeded')).all()
    h5py_compare_files_statistically(h5py_dataset_iterator,
                                        assert_allclose,
                                        tmp_path/'normal.hdf5',
                                        tmp_path/'restart.hdf5',
                                        derivative)
