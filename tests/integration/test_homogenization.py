# SPDX-License-Identifier: AGPL-3.0-or-later
import copy

import pytest
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'homogenization'


results = {}
rotation = damask.Rotation.from_random()

@pytest.mark.parametrize('homogenization',['none','SX','Taylor2','RGC'])
@pytest.mark.parametrize('phase',['Elastic','Isotropic','Phenopowerlaw'])
def test_homogenization_equivalent(res_path,tmp_path,copy_files,phase,homogenization):
    global results
    copy_files(res_path,tmp_path)

    load = 'compressionX'
    grid = 'test'
    material = 'material'
    job = f'{grid}_{load}'

    material_config = damask.ConfigMaterial.load(tmp_path/f'{material}.yaml')

    N_constituents = material_config['homogenization'][homogenization]['N_constituents']

    constituent = {'phase':phase,
                   'v':1/N_constituents,
                   'O':rotation.as_quaternion().tolist()}

    # avoid flawed error for wrong N_constituents
    for h in list(material_config['homogenization'].keys()):
        if h != homogenization: del material_config['homogenization'][h]

    material_list = {'homogenization':homogenization,
                    'constituents':[copy.deepcopy(constituent) for c in range(N_constituents)]}

    material_config['material'] = [material_list]

    material_config.save(tmp_path/f'{material}.yaml')

    stdout,stderr = damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',
                                    wd=tmp_path)
    print(f'{stdout}\n{stderr}')

    result = damask.Result(tmp_path/f'{job}.hdf5').view(increments=50).view(phases=phase)
    F_11 = np.average(result.get('F'),axis=0)[0,0]
    P_11 = np.average(result.get('P'),axis=0)[0,0]

    if phase not in results:
        results[phase]={'F':F_11,'P':P_11}

    assert 0.99 < results[phase]['F']/F_11 < 1.01 and \
           0.99 < results[phase]['P']/P_11 < 1.01


def test_homogenization_many(tmp_path,res_path,copy_files,np_rng):
    grid = damask.GeomGrid(np.zeros([2,1,1],int),np.array([2,1,1])*1e-4)
    grid.save(tmp_path/'small')

    N = np_rng.integers(100,2000)
    O = damask.Rotation.from_random(N,rng_seed=np_rng)

    m = damask.ConfigMaterial()
    m['homogenization']={'Taylor':{'N_constituents':N,'mechanical':{'type':'isostrain'}}}
    m['phase']['Elastic'] = damask.YAML.load(res_path/'material.yaml')['phase']['Elastic']
    m.material_add(homogenization='Taylor',O=O.reshape((1,-1)),phase='Elastic').save(tmp_path/'material.yaml')

    copy_files(res_path,tmp_path,['none.yaml'])
    damask.util.run(f'damask_grid -l none.yaml -g small.vti -m material.yaml --wd {tmp_path}')


def test_homogenization_fraction(tmp_path,assert_allclose,np_rng):
    grid = 'small'
    load = 'mixed'
    material = 'material'
    job = f'{grid}_{load}'

    g = damask.GeomGrid(np.zeros([2,1,1],int),np.array([2,1,1])*1e-4)
    g.initial_conditions['T'] = 150 + np_rng.random()*600
    g.save(tmp_path/f'{grid}.vti')

    v = 0.1 + np_rng.random()*0.8

    m = damask.ConfigMaterial()
    m['homogenization']={'Taylor':{'N_constituents':2,
                                   'mechanical':{'type':'isostrain','output':['F','P']},
                                   'thermal':{'type':'isotemperature','output':['T']}}}
    m['phase']['Al'] = {'lattice':'cF', 'rho':2700.0,
                        'mechanical': {'elastic': {'type':'Hooke',
                                                   'C_11':106.9e9, 'C_12':60.55e9, 'C_44':28.37e9},
                                       'output': ['F','P']},
                        'thermal':{'K_11':2.380e2,'C_p':910.0,'output':['T']}}
    m['phase']['Fe'] = {'lattice':'cI', 'rho':7874.0,
                        'mechanical': {'elastic': {'type':'Hooke',
                                                   'C_11':232.1e9, 'C_12':135.9e9, 'C_44':117.0e+9},
                                       'output': ['F','P']},
                        'thermal':{'K_11':8.055e1,'C_p':450.0,'output':['T'],
                                   'source': [{'type': 'externalheat', 'f':[1e7, 5e7], 't':[0, 10]}]}}

    m = m.material_add(homogenization='Taylor',
                       O=damask.Rotation().broadcast_to((1,2)),
                       phase=np.array(['Al','Fe']).reshape(1,2),
                       v = np.array([v,1-v]).reshape(1,2))
    m.save(tmp_path/f'{material}.yaml')

    load_case = damask.YAML(solver={'mechanical':'spectral_basic',
                                    'thermal':'spectral'},loadstep=[])

    F = [['x',0,0], [0,1.01+np_rng.random()*0.01,0.01+np_rng.random()*0.01], [0,0,'x']]
    P = [['0','x','x'], ['x','x','x'], ['x','x',0]]
    loadstep = {'boundary_conditions':{'mechanical':{'P':P,'F':F}},
                'discretization':{'t':10.,'N':20}}
    load_case['loadstep'].append(loadstep)

    load_case.save(tmp_path/f'{load}.yaml')

    damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml --wd {tmp_path} -j {job}')

    result = damask.Result(tmp_path/f'{job}.hdf5').view(increments=-1)

    T_homog = result.view(phases=False).get('T')
    P_homog = result.view(phases=False).get('P')
    F_homog = result.view(phases=False).get('F')

    T_phases = result.view(homogenizations=False).get('T')
    P_phases = result.view(homogenizations=False).get('P')
    F_phases = result.view(homogenizations=False).get('F')

    assert_allclose(np.average(F_phases['Al'],0),np.average(F_homog,0))
    assert_allclose(np.average(F_phases['Fe'],0),np.average(F_homog,0))
    assert_allclose(np.average(T_phases['Al'],0),np.average(T_homog,0))
    assert_allclose(np.average(T_phases['Fe'],0),np.average(T_homog,0))

    assert_allclose(np.average(P_phases['Al']*v + P_phases['Fe']*(1-v),0),
                    np.average(P_homog,0),atol=1,rtol=0)
