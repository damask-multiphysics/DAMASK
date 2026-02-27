# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'orientation_consistency'


def test_orientation_consistency(tmp_path,res_path,copy_files,assert_allclose):

    copy_files(res_path,tmp_path,['none.yaml'])

    data = damask.Table.load(res_path/'3D_EBSD.txt')
    data = data.set('O',data.get('O')/np.linalg.norm(data.get('O'),axis=1,keepdims=True))
    material_config = damask.ConfigMaterial.from_table(data,**{'O':'O'})
    for i,m in enumerate(material_config['material']):
        material_config['material'][i]['homogenization']='SX'
        material_config['material'][i]['constituents'][0]['phase']='tI'

    material_config['homogenization'] = {'SX': {'N_constituents': 1,
                                                'mechanical': {'type': 'pass'}}}
    C = {'C_11': 106.75e9, 'C_33': 106.75e9, 'C_12': 60.41e9,
         'C_13':  60.41e9, 'C_44':  28.34e9, 'C_66': 28.34e9}
    material_config['phase'] = {'tI': {'lattice': 'tI','c/a':(8./3.)**.5,
                                       'mechanical': {'elastic': {**{'type': 'Hooke',**C}},
                                                      'output': ['O','F_p']}}}
    material_config.save(tmp_path/'material.yaml')

    damask.GeomGrid.from_table(data,'pos','O').save(tmp_path/'3D_EBSD.vti')

    damask.util.run('damask_grid -l none.yaml -g 3D_EBSD.vti -m material.yaml -j 3D_EBSD_none',
                    wd=tmp_path)

    data = damask.Table.load(res_path/'3D_EBSD.txt')
    O_ref = data.get('O')/np.linalg.norm(data.get('O'),axis=1,keepdims=True)
    r = damask.Result(tmp_path/'3D_EBSD_none.hdf5')
    r.add_rotation('F_p')
    r.view(increments=0)
    O_cur1 = r.view(increments=0).place('O')
    O_cur2 = damask.Rotation.from_matrix(r.view(increments=0).place('R(F_p)')).as_quaternion()

    assert_allclose(O_ref,O_cur1)
    assert_allclose(O_ref,O_cur2)
