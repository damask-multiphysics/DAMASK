# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'grid_parse_loadcase'


def final(delta_t,start,end):
    return rate_eng(delta_t,start,(end-start)/sum(delta_t))

def rate_eng(delta_t,start,rate):
    v = [start]
    for d in delta_t:
        v.append(v[-1]+rate*d)
    return v[1:]

def rate_true(delta_t,start,rate):
    v = [start]
    for d in delta_t:
        v.append(v[-1]+np.dot(rate,v[-1])*d)
    return v[1:]

def time_inc(discretization):
    if 'r' not in discretization or np.isclose(discretization['r'],1.0):
        return [discretization['t']/discretization['N'] for i in range(discretization['N'])]
    else:
        return [discretization['t'] * (discretization['r']**(i)-discretization['r']**(i+1)) \
                / (1.0-discretization['r']**discretization['N']) for i in range(discretization['N'])]

def to_masked_array(lst):
    return np.array([[0.0 if i == 'x' else float(i) for i in sublst] for sublst in lst], dtype='float')

def get_target_increments(loadcase):
    count = 0
    increments = [count]
    for step in loadcase:
        N = step['discretization']['N']
        f = step.get('f_out',1)
        increments += [] if f=='none' else range(count+f,count+N+1,f)
        count += N
    return increments

def calculate_target(loadcase):
    M = [np.zeros((3,3),dtype='bool')]
    P = [np.zeros((3,3))]
    F = [np.eye(3)]

    for l in loadcase:
        delta_t = time_inc(l['discretization'])

        if   'P' in l['boundary_conditions']['mechanical']:
            P_new = final(delta_t,P[-1],to_masked_array(l['boundary_conditions']['mechanical']['P']))
            M_new = [np.array(l['boundary_conditions']['mechanical']['P']).reshape(3,3) != 'x' for _ in delta_t]
        elif 'dot_P' in l['boundary_conditions']['mechanical']:
            P_new = rate_eng(delta_t,P[-1],to_masked_array(l['boundary_conditions']['mechanical']['dot_P']))
            M_new = [np.array(l['boundary_conditions']['mechanical']['dot_P']).reshape(3,3) != 'x' for _ in delta_t]
        else:
            P_new = [np.zeros((3,3)) for _ in delta_t]
            M_new = [np.zeros((3,3),dtype='bool') for _ in delta_t]

        if 'L' in l['boundary_conditions']['mechanical']:
            F_new = rate_true(delta_t,F[-1],to_masked_array(l['boundary_conditions']['mechanical']['L']))
        elif 'dot_F' in l['boundary_conditions']['mechanical']:
            F_new = rate_eng(delta_t,F[-1],to_masked_array(l['boundary_conditions']['mechanical']['dot_F']))
        elif 'F' in l['boundary_conditions']['mechanical']:
            F_new = final(delta_t,F[-1],to_masked_array(l['boundary_conditions']['mechanical']['F']))

        P += P_new
        M += M_new
        F += F_new

    incs = get_target_increments(loadcase)
    return([F[i] for i in incs],[P[i] for i in incs],[M[i] for i in incs])


s = (1+np.random.rand())*1e-3
loadcases = [
             [{'boundary_conditions': {'mechanical': {'dot_F': [[ 0,   s,   0],  [  0,  0,  0],  ['x','x','x']],
                                                      'P':     [['x', 'x', 'x'], [ 'x','x','x'], [ 0,  0,  0 ]]}},
               'discretization': {'t': 100,'N': 100},'f_out': np.random.randint(2,50)}],
             [{'boundary_conditions': {'mechanical': {'dot_F': [['x',  0,   0],  [  s, 'x', 0],  [  0,   0, 'x']],
                                                      'P':     [[ 0,  'x', 'x'], [ 'x', 0, 'x'], [ 'x', 'x', 0 ]]}},
               'discretization': {'t': 100,'N': 100,'r': 1.0},'f_out': 10}],
             [{'boundary_conditions': {'mechanical': {'dot_F': [['x',  0,   0],  [  0,   0,   s], [ 0,   0, 'x']],
                                                      'P':     [[ 0,  'x', 'x'], [ 'x', 'x', 'x'],['x', 'x', 0 ]]}},
               'discretization': {'t': 20, 'N': 20, 'r': 1.1}}],
             [{'boundary_conditions': {'mechanical': {'P':     [[ 0,   0,   0],  [ 'x', 'x', 'x'], [ 'x', 'x', 'x']],
                                                      'L':     [['x', 'x', 'x'], [  0,   0,   0],  [ 1e-3, 0 ,   s]]}},
               'discretization': {'t': 20, 'N': 20, 'r': 0.95}}],
             [{'boundary_conditions': {'mechanical': {'dot_F': [['x', 0, 'x'],  [ 0,-1e-3,2e-3], [ 0,  0,  0]],
                                                      'P':     [[ 0, 'x', 0],   ['x','x', 'x'],  ['x','x','x']]}},
               'discretization': {'t': 30,'N': 150},'f_out': 10}],
             [{'boundary_conditions': {'mechanical': {'F':     [[1.1,0,0],
                                                                [0, float(np.sqrt(1/1.1)),0],
                                                                [ 0, 0,float(np.sqrt(1/1.1))]]}},
               'discretization': {'t': 50,'N': 30}}],
             [{'boundary_conditions': {'mechanical': {'dot_F': [['x',  0,   0],   [ 0,   0,   0],   [ 0,   0, 'x']],
                                                      'P':     [[4e7, 'x', 'x'],  ['x', 'x', 'x'],  ['x', 'x', 0 ]]}},
               'discretization': {'t': 40,'N': 200}},
              {'boundary_conditions': {'mechanical': {'dot_F': [['x',  0,   0],   [1e-2,   0,   0], [   0,   0, 'x']],
                                                      'P':     [[1e2, 'x', 'x'],  [ 'x', 'x',  'x'],[  'x', 'x', 0 ]]}},
               'discretization': {'t': 40,'N': 200}}],
             [{'boundary_conditions': {'mechanical': {'dot_F': [['x',  0,   0],   [ 0,   0,  -s],  [ 0,   0, 'x']],
                                                      'P':     [[1e7, 'x', 'x'],  ['x', 'x', 'x'], ['x', 'x', 0 ]]}},
               'discretization': {'t': 20,'N': 20}},
              {'boundary_conditions': {'mechanical': {'dot_F': [['x',  0,   0],   [ 0,   0,  -s],  [ 0,   0, 'x']],
                                                      'P':     [[1e7, 'x', 'x'],  ['x', 'x', 'x'], ['x', 'x', 0 ]]}},
               'discretization': {'t': 60,'N': 20}}],
             [{'boundary_conditions': {'mechanical': {'dot_F': [['x', 0, 'x'],  [ 0,-1e-3,2e-3], [ 0,  0,  0]],
                                                      'dot_P': [[ 0, 'x', 0],   ['x','x', 'x'],  ['x','x','x']]}},
               'discretization': {'t': 30,'N': 150},'f_out': np.random.randint(4,40)}],
             [{'boundary_conditions': {'mechanical': {'dot_F': [['x',  0,   0],   [ 0,   0,  -s],  [ 0,   0, 'x']],
                                                      'P':     [[1e7, 'x', 'x'],  ['x', 'x', 'x'], ['x', 'x', 0 ]]}},
               'discretization': {'t': 20,'N': 60}, 'f_out': 8},
              {'boundary_conditions': {'mechanical': {'dot_F': [['x',  0,   0],   [ 0,   0,  -s],  [ 0,   0, 'x']],
                                                      'P':     [[1e7, 'x', 'x'],  ['x', 'x', 'x'], ['x', 'x', 0 ]]}},
               'discretization': {'t': 20,'N': 20}, 'f_out': 'none'},
              {'boundary_conditions': {'mechanical': {'dot_F': [['x',  0,   0],   [ 0,   0,  -s],  [ 0,   0, 'x']],
                                                      'P':     [[1e7, 'x', 'x'],  ['x', 'x', 'x'], ['x', 'x', 0 ]]}},
               'discretization': {'t': 60,'N': 22}},
              {'boundary_conditions': {'mechanical': {'dot_F': [['x',  0,   0],   [ 0,   0,  -s],  [ 0,   0, 'x']],
                                                      'P':     [[1e7, 'x', 'x'],  ['x', 'x', 'x'], ['x', 'x', 0 ]]}},
               'discretization': {'t': 60,'N': 22}, 'f_out': 30},
              {'boundary_conditions': {'mechanical': {'dot_F': [['x',  0,   0],   [ 0,   0,  -s],  [ 0,   0, 'x']],
                                                      'P':     [[1e7, 'x', 'x'],  ['x', 'x', 'x'], ['x', 'x', 0 ]]}},
               'discretization': {'t': 60,'N': 20}, 'f_out': 'none'},
              {'boundary_conditions': {'mechanical': {'dot_F': [['x',  0,   0],  [  s, 'x', 0],  [  0,   0, 'x']],
                                                      'P':     [[ 0,  'x', 'x'], [ 'x', 0, 'x'], [ 'x', 'x', 0 ]]}},
               'discretization': {'t': 100,'N': 100,'r': 1.1},'f_out': 100}],
            ]

@pytest.mark.parametrize('solver',['spectral_basic',
                                   pytest.param('FEM',marks=pytest.mark.xfail),
                                   'spectral_polarization',
                                   'spectral_Galerkin'])
@pytest.mark.parametrize('loadcase',loadcases)
def test_grid_parse_loadcase(tmp_path,res_path,copy_files,
                             solver,loadcase):
    load = 'load'
    grid = 'simple'
    job = f'{grid}_{load}'

    copy_files(res_path,tmp_path)

    load_config = damask.YAML({'solver': {'mechanical': solver}, 'loadstep': loadcase})
    load_config.save(tmp_path/f'{load}.yaml')

    cmd = f'damask_grid -l {load}.yaml -g {grid}.vti -m material.yaml -j {job}'
    if solver == 'FEM':
        petsc_options = '-pc_type none -snes_type newtontr'
        damask.YAML({'solver': {'grid': {'mechanical': {'PETSc_options':petsc_options}}}}).save(tmp_path/'numerics.yaml')
        cmd += ' -n numerics.yaml'

    out, _ = damask.util.run(cmd,wd=tmp_path)
    assert solver in out

    F_target,P_target,mask = calculate_target(loadcase)
    r = damask.Result(tmp_path/f'{job}.hdf5')

    for i,F in enumerate([np.average(F_field,axis=0) for F_field in r.get('F').values()]):
        F_cur = np.ma.masked_array(F,mask[i])
        F_ref = np.ma.masked_array(F_target[i],mask[i])
        assert np.ma.allclose(F_cur,F_ref,atol=1.e-2,rtol=0)

    for i,P in enumerate([np.average(P_field,axis=0) for P_field in r.get('P').values()]):
        P_cur = np.ma.masked_array(P,np.logical_not(mask[i]))
        P_ref = np.ma.masked_array(P_target[i],np.logical_not(mask[i]))
        assert np.all(np.abs(P_cur-P_ref)<max(np.max(np.abs(P))*1e-3,1e3)) in [True,np.ma.masked]
