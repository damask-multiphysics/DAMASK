# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import h5py

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'mesh_element_lib'


@pytest.mark.parametrize('dimension',['3D','2D'])
@pytest.mark.parametrize('p_s,p_i',
                         [(1,1),
                          (2,2),(2,1),
                          (3,3),(3,2),
                          (4,4),(4,3),
                          (5,5),(5,4)])
@pytest.mark.xfail
def test_mesh_element_lib(res_path,tmp_path,copy_files,update,assert_allclose,dimension,p_i,p_s):

    load = 'tensionX'
    mesh = f'bicrystal_{dimension}'
    material = 'material'
    job = f'{mesh}_{load}'
    copy_files(res_path,tmp_path,[f'{mesh}.msh',f'{material}.yaml','tensionX.yaml'])

    numerics_config = damask.YAML.load(res_path/'numerics.yaml')
    numerics_config['solver']['mesh']['p_i'] = p_i
    numerics_config['solver']['mesh']['p_s'] = p_s
    numerics_config.save(tmp_path/'numerics.yaml')

    damask.util.run(f'DAMASK_mesh -l {load}.yaml -g {mesh}.msh -m {material}.yaml -j {job} -n numerics.yaml',wd=tmp_path)

    with h5py.File(res_path/'Reference.hdf5','a') as ref, h5py.File(tmp_path/f'{job}.hdf5','r') as cur:
        for root, datasets in zip(('/','increment_200'),(('x_n','x_p'),('u_n','u_p'))):
            for dataset in datasets:
                loc_cur = f'{root}/geometry/{dataset}'
                loc_ref = f'/{dimension}-p_s{p_s}-p_i{p_i}/{dataset}'
                if update:
                    if loc_ref not in ref: ref.create_dataset(loc_ref, data=cur[loc_cur])
                    ref[loc_ref][()]=cur[loc_cur][()]
                assert_allclose(ref[loc_ref][()],cur[loc_cur][()],rtol=5e-4,atol=5e-6)
