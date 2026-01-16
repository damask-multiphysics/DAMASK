# SPDX-License-Identifier: AGPL-3.0-or-later
import shutil
import tempfile
from pathlib import Path

import pytest
import numpy as np
import h5py

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'state_integrators'

ref_file = None


@pytest.mark.parametrize('increments',[150,50,40])
@pytest.mark.parametrize('integrator',['FPI','Euler','Euler_adaptive','RK4','RKCK45'])
def test_state_integrators_detect_changes(integrator,increments,
                                          res_path,tmp_path,copy_files,h5py_dataset_iterator,assert_allclose):
    global ref_file
    load='tensionX'
    grid='test'
    material = 'material'
    job=f'{grid}_{load}'

    copy_files(res_path,tmp_path,[f'{grid}.vti',f'{material}.yaml'])

    numerics_config = damask.YAML.load(res_path/'numerics.yaml')
    numerics_config['phase']['mechanical']['plastic']['integrator_state'] = integrator
    numerics_config.save(tmp_path/'numerics.yaml')

    l = damask.LoadcaseGrid.load(res_path/f'{load}.yaml')
    l['loadstep'][1]['discretization']['N']=increments
    l['loadstep'][1]['f_out']=increments
    l.save(tmp_path/f'{load}.yaml')

    damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml '+
                    f'-n numerics.yaml -j {job}',wd=tmp_path)

    if ref_file is None:
        print('copied reference')
        ref_file = Path(tempfile.mkdtemp())/f'{job}.hdf5'
        shutil.copyfile(tmp_path/f'{job}.hdf5',ref_file)

    with h5py.File(ref_file,'r') as ref, \
            h5py.File(tmp_path/f'{job}.hdf5','r+') as cur:
        if increments != 150: cur.move(f'increment_{increments+10}','increment_160')
        for (path,dset) in h5py_dataset_iterator(ref):
            if len(dset.dtype) > 0:                                                             # cell_to
                for e in dset.dtype.fields:
                    assert np.all(dset[e]==cur[path][e])
            elif dset.dtype.char == 'S':                                                        # setup
                pass
            else:
                atol=max(5e-3*max(np.amax(np.abs(dset)),np.amax(np.abs(cur[path]))),1e-12)
                assert_allclose(dset,cur[path],rtol=.0,atol=atol)
