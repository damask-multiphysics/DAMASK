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
    return res_path_base/'grid_MPI'


@pytest.mark.parametrize('grid',['3D','2D'])
@pytest.mark.parametrize('n',[2,3])
@pytest.mark.parametrize('solver',['spectral_basic','spectral_polarization','FEM','spectral_Galerkin'])
def test_grid_MPI(res_path,tmp_path,copy_files,h5py_dataset_iterator,mpi_launcher,assert_allclose,
                  solver,n,grid):

    load = 'tensionY'
    material = 'material'
    job = f'{grid}_{load}'

    for mode in ['serial','parallel']:
        cwd = tmp_path/mode

        os.mkdir(cwd)
        copy_files(res_path,cwd)
        load_config = damask.YAML.load(cwd/f'{load}.yaml')
        load_config['solver']['mechanical'] = solver
        load_config.save(cwd/f'{load}.yaml')
        out,_ = damask.util.run((f'{mpi_launcher} -n {n} ' if mode == 'parallel' else '') +
                                 f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',
                                wd=cwd,env=env)
        assert solver in out

    serial = h5py.File(tmp_path/f'serial/{job}.hdf5','r')
    parallel = h5py.File(tmp_path/f'parallel/{job}.hdf5','r')
    for (path,dset) in h5py_dataset_iterator(serial):
        if len(dset.dtype) > 0:                                                                     # cell_to
            for e in dset.dtype.fields:
                assert np.all(dset[e]==parallel[path][e])
        elif dset.dtype.char == 'S':                                                                # setup
            dset == parallel[path]
        else:
            atol=max(5e-2*max(np.amax(np.abs(dset)),np.amax(np.abs(parallel[path]))),1e-12)
            assert_allclose(dset,parallel[path],rtol=.0,atol=atol)
