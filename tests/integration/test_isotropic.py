# SPDX-License-Identifier: AGPL-3.0-or-later
import os

import pytest
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'isotropic'


@pytest.fixture
def mat_configs(np_rng):
    """Random material configuration (isotropic/phenopowerlaw)."""
    C_11 = (15+np_rng.random()*10)*1e10
    C_12 = C_11 - (5+np_rng.random()*5)*1e10
    C_44 = (C_11-C_12)*.5 + (np_rng.random()-0.5)*2e10

    E = C_44*(3*C_12+2*C_44)/(C_44+C_12)                                                            # isotropic stiffness

    a = 1 + np_rng.random()*4
    n = 5 + np_rng.random()*30
    dot_gamma_0 = (1e5**np_rng.random())*1e-4
    xi_0 = E*.002 + (1 + np_rng.random()*10)*1e7
    xi_inf = xi_0*(1+np_rng.random()*5)
    h_0 = (1 + np_rng.random()*200)*1e7
    h = float(np_rng.choice([1.0,2.5,3.0]))
    h_sl_sl = {1.0:([1]+[0]*6), 2.5:([1]*7), 3.0:([1,1]+[1.4]*5)}[h]

    homog = {'SX': {'N_constituents':1, 'mechanical':{'type':'pass'}}}

    elastic = {'type':'Hooke', 'C_11':C_11, 'C_12':C_12, 'C_44':C_44}
    mech = {'lattice': 'cF', 'output':['F','P','F_p','F_e','L_p'], 'elastic': elastic}

    isotropic = {'type':'isotropic', 'a':a, 'n':n, 'dot_gamma_0':dot_gamma_0,
                 'xi_0':xi_0, 'xi_inf':xi_inf, 'h_0':h_0, 'M':2.5, 'h':h,
                 'output':['xi']}
    pheno = {'type':'phenopowerlaw', 'a_sl':[a], 'n_sl':[n], 'dot_gamma_0_sl':[dot_gamma_0],
             'xi_0_sl':[xi_0], 'xi_inf_sl':[xi_inf], 'h_0_sl-sl':[h_0],
             'N_sl':[12], 'h_sl-sl':h_sl_sl}

    mat = damask.ConfigMaterial()
    mat['homogenization'] = homog

    configs = dict()
    for model,config in zip(['isotropic','pheno'],[isotropic,pheno]):
        configs[model] = mat.copy()
        configs[model]['phase'][model] = {'lattice':'cF', 'mechanical':mech.copy()}
        configs[model]['phase'][model]['mechanical']['plastic'] = config.copy()

    return configs


@pytest.mark.parametrize('load',['tensionZ','plane_strain_compression'])
def test_phenopowerlaw_equivalence(res_path,tmp_path,copy_files,mat_configs,np_rng,load):

    cell = 11
    grid = damask.GeomGrid(np.arange(cell**3).reshape([cell]*3),np.ones(3)*1e-5)
    rot = damask.Rotation.from_random(cell**3,rng_seed=np_rng)

    stress = []
    for model,mat_config in mat_configs.items():
        os.mkdir(tmp_path/model)

        grid.save(tmp_path/model/f'{cell**3}.vti')
        mat_config = mat_config.material_add(phase=model,homogenization='SX',O=rot)
        mat_config.save(tmp_path/model/'material.yaml')
        copy_files(res_path,tmp_path/model)

        damask.util.run(f'damask_grid -l {load}.yaml -g {cell**3}.vti -m material.yaml '\
                        f'--wd {tmp_path/model} --job {cell**3}_{load}')

        r = damask.Result(tmp_path/model/f'{cell**3}_{load}.hdf5')
        stress.append(np.array([np.average(P[:,2,2]) for P in r.get('P').values()]))

    assert np.average(np.abs((stress[0]-stress[1])/(stress[0]+stress[1])*2)) < .1


def test_analytic_reference(res_path,tmp_path,copy_files,mat_configs,np_rng):
    mat_config = mat_configs['isotropic']
    C_11 = mat_config['phase']['isotropic']['mechanical']['elastic']['C_11']
    C_12 = mat_config['phase']['isotropic']['mechanical']['elastic']['C_12']
    C_44 = (C_11-C_12)*.5
    mat_config['phase']['isotropic']['mechanical']['elastic']['C_44'] = C_44                        # isotropic, simplifies calculation of E
    M = 1+np_rng.random()*3
    mat_config['phase']['isotropic']['mechanical']['plastic']['M'] = M

    E = C_44*(3*C_12+2*C_44)/(C_44+C_12)

    a = mat_config['phase']['isotropic']['mechanical']['plastic']['a']
    n = mat_config['phase']['isotropic']['mechanical']['plastic']['n']
    dot_gamma_0 = mat_config['phase']['isotropic']['mechanical']['plastic']['dot_gamma_0']

    xi_0 = mat_config['phase']['isotropic']['mechanical']['plastic']['xi_0']
    xi_inf = mat_config['phase']['isotropic']['mechanical']['plastic']['xi_inf']

    h_0 = mat_config['phase']['isotropic']['mechanical']['plastic']['h_0']
    h = mat_config['phase']['isotropic']['mechanical']['plastic']['h']


    grid = damask.GeomGrid(np.zeros([2,1,1],int),np.ones(3)*1e-5)
    grid.save(tmp_path/'homogeneous.vti')

    mat_config = mat_config.material_add(phase='isotropic',homogenization='SX',O=damask.Rotation())
    mat_config.save(tmp_path/'material.yaml')

    copy_files(res_path,tmp_path)

    damask.util.run(f'damask_grid -l tensionX.yaml -g homogeneous.vti -m material.yaml --wd {tmp_path} '\
                        '--job homogeneous_tensionX')

    r = damask.Result(tmp_path/'homogeneous_tensionX.hdf5')
    sim = {}
    for out in ['F','P','F_e','F_p','L_p']:
        sim[out] = np.array([np.average(o[:,0,0]) for o in r.get(out).values()])
    sim['xi'] = np.array([np.average(o) for o in r.get('xi').values()])

    load = damask.YAML.load(tmp_path/'tensionX.yaml')
    dot_F = load['loadstep'][0]['boundary_conditions']['mechanical']['dot_F'][0][0]
    t = load['loadstep'][0]['discretization']['t']
    N = load['loadstep'][0]['discretization']['N']
    Delta_t = t/N

    F = [1]
    F_p = [1]
    xi = [xi_0]
    S = [0]

    dot_xi = 0
    dot_F_p = 0

    for inc in range(N):

        F.append(F[-1] + dot_F*Delta_t)

        F_p_guess = F_p[-1] + dot_F_p*Delta_t
        xi_guess = xi[-1] + dot_xi*Delta_t

        for i in range(5):                                                                          # attempt of FPI

            F_e_guess = F[-1]/F_p_guess
            S_guess = E * (F_e_guess**2-1)*.5                                                       # 1D Hooke
            norm_Mp_dev = (S_guess/1.5**.5)                                                         # S = M_p for F_i = I_3

            dot_gamma = dot_gamma_0 * (1.5**.5 * norm_Mp_dev/(M*xi_guess))**n
            L_p_guess = dot_gamma/1.5**.5
            dot_F_p = L_p_guess*F_p[-1]
            dot_xi = dot_gamma * h_0 * np.abs(1-(xi_guess/xi_inf))**a * np.sign(1.-xi_guess/xi_inf) *h

            xi_guess = xi[-1] + dot_xi*Delta_t
            F_p_guess = F_p[-1] + dot_F_p*Delta_t

        S.append(S_guess)
        F_p.append(F_p_guess)
        xi.append(xi_guess)

    F_e = np.array(F)/np.array(F_p)
    P = np.array(S)/np.array(F_p)**2*np.array(F)

    assert np.average(np.abs((sim['F']-np.array(F))/sim['F'])) < 1e-12
    assert np.average(np.abs((sim['F_e']-np.array(F_e))/sim['F_e'])) < 5e-4
    assert np.average(np.abs((sim['F_p']-np.array(F_p))/sim['F_p'])) < 5e-4
    assert np.average(np.abs((sim['P'][1:]-np.array(P[1:]))/sim['P'][1:])) < 5e-2
    assert np.average(np.abs((sim['xi']-np.array(xi))/sim['xi'])) < 5e-2
