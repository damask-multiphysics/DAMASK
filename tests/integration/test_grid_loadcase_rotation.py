# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np
import h5py

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'grid_loadcase_rotation'


equivalent_loadcases = [{'N':250,'freq':25,'t':500,'deformation_type':'dot_F',
                         'deformation':[[1e-4,0,0,  0, 'x',0,  0,0,0],
                                        [ 'x',0,0,  0,1e-4,0,  0,0,0],
                                        [ 'x',0,0,  0,1e-4,0,  0,0,0],
                                        [   0,0,0,  0, 'x',0,  0,0,1e-4]],
                         'rotation':[[1,0,0,0],[0,0,1,90],[0,0,-1,90],[0,1,0,90]]
                         },
                        {'N':250,'freq':25,'t':500,'deformation_type':'L',
                         'deformation':[[   0,0,0,    0,1e-4,0,  'x','x','x'],
                                        [1e-4,0,0,    0,   0,0,  'x','x','x'],
                                        [1e-4,0,0,    0,   0,0,  'x','x','x'],
                                        [   0,0,0,  'x', 'x','x',  0,  0, 1e-4]],
                         'rotation':[[1,0,0,0],[0,0,1,90],[0,0,-1,90],[1,0,0,90]]
                         },
                        {'N':250,'freq':25,'t':500,'deformation_type':'F',
                         'deformation':[[  'x',0,0, 0,'x',0,  0,0,1.2],
                                        [  1.2,0,0, 0,'x',0,  0,0,'x'],
                                        [  1.2,0,0, 0,'x',0,  0,0,'x'],
                                        [  'x',0,0, 0,1.2,0,  0,0,'x'],
                                        [  'x',0,0, 0,1.2,0,  0,0,'x']],
                         'rotation':[[1,0,0,0],[0,1,0,90],[0,-1,0,90],[1,0,0,90],[-1,0,0,90]]
                         }
                        ]

def loadcase_string(loadcases):
    load_case = []
    for d,R in zip(loadcases['deformation'],loadcases['rotation']):
        P = [0 if m == 'x' else 'x' for m in d]
        load_case.append([{'boundary_conditions': {'mechanical': {loadcases['deformation_type']: [d[0:3],d[3:6],d[6:9]],
                                                                  'P': [P[0:3],P[3:6],P[6:9]],
                                                                  'R': R}},
                           'discretization': {'t': loadcases['t'], 'N': loadcases['N']},
                           'f_out': loadcases['freq']}])
    return load_case

@pytest.mark.parametrize('loadcases',equivalent_loadcases)
@pytest.mark.parametrize('solver',['spectral_basic','spectral_polarization','FEM'])
def test_grid_loadcase_rotation(res_path,tmp_path,copy_files,h5py_dataset_iterator,assert_allclose,
                                loadcases,solver):

    grid = 'simple'
    material = 'material'
    cmd_base = f'damask_grid -g {grid}.vti -m {material}.yaml'

    copy_files(res_path,tmp_path)
    if solver == 'FEM':
        petsc_options = '-pc_type none -snes_type newtontr'
        damask.YAML({'solver': {'grid': {'mechanical': {'PETSc_options':petsc_options}}}}).save(tmp_path/'numerics.yaml')
        cmd_base += ' -n numerics.yaml'

    for i,s in enumerate(loadcase_string(loadcases)):
        load_config = damask.YAML({'solver': {'mechanical': solver},'loadstep': s})
        load_config.save(tmp_path/f'{i}.yaml')

        job = f'{grid}_{i}'
        out, _ = damask.util.run(cmd_base+f' -j {job} -l {i}.yaml',wd=tmp_path)
        assert solver in out

        with h5py.File(tmp_path/f'{grid}_0.hdf5','r') as ref, \
             h5py.File(tmp_path/f'{job}.hdf5','r') as cur:
            for (path,dset) in h5py_dataset_iterator(ref):
                if len(dset.dtype) > 0:                                                             # cell_to
                    for e in dset.dtype.fields:
                        assert np.all(dset[e]==cur[path][e])
                elif dset.dtype.char == 'S':                                                        # setup
                    pass
                else:
                    if 'geometry/u' not in path and solver == 'spectral_polarization': # displacement broken
                        atol=max(1e-8*max(np.amax(np.abs(dset)),np.amax(np.abs(cur[path]))),1e-12)
                        assert_allclose(dset,cur[path],rtol=.0,atol=atol)
