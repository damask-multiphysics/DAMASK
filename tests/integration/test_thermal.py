# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'thermal'


@pytest.mark.parametrize('lattice',['cI','cF','hP','tI'])
@pytest.mark.parametrize('case,T_0,result',[('adiabatic',         400, [400,900]),
                                            ('fast_convection',   340, [356.6015625,356.6015625])])
def test_conduction(res_path,tmp_path,copy_files,assert_allclose,
                    lattice,case,T_0,result):
    grid = 'inclusion'
    load = 'no_deformation'
    material = 'material'
    job  = f'{grid}_{load}'
    copy_files(res_path,tmp_path,[f'{load}.yaml',f'{material}.yaml'])

    mat = damask.ConfigMaterial.load(tmp_path/f'{material}.yaml')
    for phase in ['heatsource', 'matrix']:
        mat['phase'][phase]['lattice'] = lattice
        if lattice in ['hP','tI']:
            mat['phase'][phase]['c/a'] = (8./3.)**.5
        if case == 'adiabatic':
            mat['phase'][phase]['thermal']['K_11'] = mat['phase'][phase]['thermal']['K_33'] = 0.
        if case == 'fast_convection':
            mat['phase'][phase]['thermal']['K_11'] = mat['phase'][phase]['thermal']['K_33'] = 1.e30


    mat.save(tmp_path/f'{material}.yaml')

    g = damask.GeomGrid.load(res_path/f'{grid}.vti')
    g.initial_conditions['T'] = T_0
    g.save(tmp_path/grid)


    damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',wd=tmp_path)
    m = damask.GeomGrid.load(tmp_path/grid).material.flatten(order='F')
    r = damask.Result(tmp_path/f'{job}.hdf5')
    T = r.view(increments=10).get('T')
    assert_allclose(T,np.where(m==0,result[0],result[1]))


@pytest.mark.parametrize('lattice',['cI','cF','hP','tI'])
def test_linear_expansion(res_path,tmp_path,copy_files,assert_allclose,np_rng,
                          lattice):
    grid = 'simple'
    load = 'no_stress'
    material = 'material'
    job  = f'{grid}_{load}'
    copy_files(res_path,tmp_path,[f'{load}.yaml',f'{material}.yaml'])

    g = damask.GeomGrid(np.zeros((2,2,2)),np.ones(3)*1e-6)
    g.initial_conditions['T'] = np_rng.integers(1000)
    g.save(tmp_path/grid)

    mat = damask.ConfigMaterial.load(tmp_path/f'{material}.yaml')
    mat['phase']['matrix']['lattice'] = lattice
    if lattice in ['hP','tI']:
        mat['phase']['matrix']['c/a'] = (8./3.)**.5
    mat['phase']['matrix']['thermal']['K_11'] = mat['phase']['matrix']['thermal']['K_33'] = 0.0
    mat['phase']['heatsource']['thermal']['K_11'] = mat['phase']['heatsource']['thermal']['K_33'] = 0.0
    mat['phase']['matrix']['mechanical']['eigen'] = [{'Alpha_11':1e-5, 'Alpha_33': 0,
                                                        'T_ref':np_rng.integers(1,1000), 'type': 'thermalexpansion'}]
    mat['phase']['matrix']['thermal']['source'] = [{'type':'externalheat', 'f':[1,1], 't':[0,100]}]

    mat['material'][0]['constituents'][0]['O'] = damask.Rotation.from_random(rng_seed=np_rng)
    mat.save(tmp_path/f'{material}.yaml')

    damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',wd=tmp_path)
    r = damask.Result(tmp_path/f'{job}.hdf5')
    r.add_strain()
    l = {}
    for v in ['min','mid','max']:
        r.add_eigenvalue('epsilon_V^0.0(F)',v)
        l[v] = np.array([np.average(x) for x in r.get(f'lambda_{v}(epsilon_V^0.0(F))').values()])

    assert_allclose(l['max'],l['mid'])
    assert_allclose(l['min'],l['mid'] if lattice.startswith('c') else 0)

    assert_allclose(l['max'][-1],1e-3,atol=5e-5,rtol=0)

def test_1D_expansion(res_path,tmp_path,copy_files,assert_allclose,np_rng):
    def integral_a(a_0,a_1,a_2,T_ref,T):
        return a_0*T + a_1/2.0*(T-T_ref)**2 + a_2/3.0*(T-T_ref)**3

    grid = 'simple'
    load = 'no_stress'
    material = 'material'
    job  = f'{grid}_{load}'
    copy_files(res_path,tmp_path,[f'{load}.yaml',f'{material}.yaml'])

    mat = damask.ConfigMaterial.load(tmp_path/f'{material}.yaml')
    mat['phase']['heatsource']['lattice'] = 'hP'
    mat['phase']['heatsource']['c/a'] = (8./3.)**.5

    a_0,a_1,a_2 = map(float,np.array([1.e-3,1.e-5,1.e-7]) \
                + (np_rng.random(3) - 0.5) * np.array([2e-4, 2e-6, 2e-8]))
    T_0 = np_rng.integers(200,500)
    T_ref = np_rng.integers(T_0,T_0+100)
    f = np_rng.random()+0.2
    mat['phase']['matrix']['thermal']['K_11'] = mat['phase']['matrix']['thermal']['K_33'] = 0.0
    mat['phase']['heatsource']['thermal']['K_11'] = mat['phase']['heatsource']['thermal']['K_33'] = 0.0
    mat['phase']['heatsource']['mechanical']['eigen'] = [{'Alpha_11':0,
                                                          'Alpha_33': a_0, 'Alpha_33,T':a_1, 'Alpha_33,T^2':a_2,
                                                          'T_ref':T_ref, 'type': 'thermalexpansion'}]
    mat['phase']['heatsource']['thermal']['source'] = [{'type':'externalheat', 'f':[f]*2, 't':[0,100]}]
    mat.save(tmp_path/f'{material}.yaml')

    g = damask.GeomGrid(np.ones((2,2,2)),np.ones(3)*1e-6)
    g.initial_conditions['T'] = T_0
    g.save(tmp_path/grid)

    damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',wd=tmp_path)
    r = damask.Result(tmp_path/f'{job}.hdf5').view(increments=-1)
    r.add_strain()
    T_1 = np.average(r.get('T'))
    eps = np.average(r.get('epsilon_V^0.0(F)')[:,2,2])

    eps_analytic = np.exp(integral_a(a_0,a_1,a_2,T_ref,T_1)-integral_a(a_0,a_1,a_2,T_ref,T_0)) -1
    assert_allclose(eps_analytic,eps,atol=1e-2,rtol=5e-3)

@pytest.mark.parametrize('lattice',['cI','cF','hP','tI'])
def test_temperature_dependent_stiffness(res_path,tmp_path,copy_files,assert_allclose,np_rng,
                                         lattice):
    load = 'tensionX'
    material = 'material'
    copy_files(res_path,tmp_path,[f'{load}.yaml',f'{material}.yaml'])

    T_ref = np_rng.integers(300,900)
    T_l = T_ref - np_rng.integers(20,200)
    T_h = T_ref + np_rng.integers(20,200)
    T = T_l if T_h%2==0 else T_h

    f_l = 0.9 - np_rng.random()*0.1
    f_h = 1.1 + np_rng.random()*0.1

    coeff = np.linalg.inv(np.array([[1,T_l,T_l**2],[1,T_ref,T_ref**2],[1,T_h,T_h**2]]))@np.array([f_l,1,f_h])

    mat = damask.ConfigMaterial.load(tmp_path/f'{material}.yaml')
    for phase in ['matrix','heatsource']:
        mat['phase'][phase]['lattice'] = lattice
        if lattice in ['hP','tI']:
            mat['phase'][phase]['c/a'] = (8./3.)**.5
        mat['phase'][phase]['thermal']['K_11'] = mat['phase'][phase]['thermal']['K_33'] = 0.0
        mat['phase'][phase]['mechanical']['elastic']['T_ref'] = T_ref

    for label in ['C_11','C_12','C_13','C_33','C_44','C_66']:
        C_xx = mat['phase']['matrix']['mechanical']['elastic'][label]
        mat['phase']['matrix']['mechanical']['elastic'][label+',T'] = C_xx*coeff[1].item()
        mat['phase']['matrix']['mechanical']['elastic'][label+',T^2'] = C_xx*coeff[2].item()
        C_xx_f = C_xx + C_xx*coeff[1].item()*(T-T_ref) + C_xx*coeff[2].item()*(T-T_ref)**2
        mat['phase']['heatsource']['mechanical']['elastic'][label] = C_xx_f

    del(mat['phase']['heatsource']['thermal']['source'])
    mat.save(tmp_path/f'{material}.yaml')

    v = damask.GeomGrid(np.zeros((2,2,2)),np.ones(3)*1e-6)
    v.initial_conditions['T'] = T
    v.save(tmp_path/'variable')

    f = damask.GeomGrid(np.ones((2,2,2)),np.ones(3)*1e-6)
    f.initial_conditions['T'] = T_ref
    f.save(tmp_path/'fixed')

    damask.util.run(f'DAMASK_grid -l {load}.yaml -g variable.vti -m {material}.yaml -j variable_{load}',
                    wd=tmp_path)
    damask.util.run(f'DAMASK_grid -l {load}.yaml -g fixed.vti -m {material}.yaml -j fixed_{load}',
                    wd=tmp_path)

    r_v = damask.Result(f'{tmp_path}/variable_{load}.hdf5').view(increments=-1)
    P_v = np.average(r_v.place('P'),0)
    T_v = np.average(r_v.place('T'),0)
    r_v = damask.Result(f'{tmp_path}/fixed_{load}.hdf5').view(increments=-1)
    P_f = np.average(r_v.place('P'),0)
    T_f = np.average(r_v.place('T'),0)
    assert_allclose(P_f,P_v,atol=1.e-3)
    assert np.isclose(T_f,T_ref) and np.isclose(T_v,T)


def test_heat_capacity(res_path,tmp_path,np_rng):

    grid = 'simple'
    load = 'no_stress'
    material = 'material'
    job  = f'{grid}_{load}'

    l = damask.LoadcaseGrid.load(res_path/f'{load}.yaml')
    t = np_rng.integers(80,150)
    l['loadstep'][0]['discretization']['t'] = t
    l['loadstep'][0]['discretization']['N'] = t
    l.save(tmp_path/f'{load}.yaml')

    mat = damask.ConfigMaterial.load(res_path/f'{material}.yaml')
    del mat['phase']['matrix']['thermal']
    del mat['phase']['heatsource']['thermal']['K_33']
    mat['phase']['heatsource']['thermal']['K_11'] = 1.0
    mat['phase']['heatsource']['lattice'] = 'cI'

    T_0 = np_rng.integers(200,500)
    T_ref = np_rng.integers(T_0-100,T_0+100)
    f = np_rng.random()+0.2
    C_p = (np_rng.random(3)+0.2)*np.array([1.,1e-2,1e-3])
    rho = np_rng.random()+0.2

    mat['phase']['heatsource']['rho'] = rho
    mat['phase']['heatsource']['thermal']['C_p'] = C_p[0]
    mat['phase']['heatsource']['thermal']['C_p,T'] = C_p[1]
    mat['phase']['heatsource']['thermal']['C_p,T^2'] = C_p[2]
    mat['phase']['heatsource']['thermal']['T_ref'] = T_ref
    mat['phase']['heatsource']['thermal']['source'] = [{'type':'externalheat', 'f':[f]*2, 't':[0,t]}]
    mat.save(tmp_path/f'{material}.yaml')

    g = damask.GeomGrid(np.ones((2,2,2)),np.ones(3)*1e-4)
    g.initial_conditions['T'] = T_0
    g.save(tmp_path/grid)

    damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',wd=tmp_path)
    T_1_sim = np.average(damask.Result(tmp_path/f'{job}.hdf5').view(increments=-1).get('T'))

    c_0 = - (f*t/rho + (C_p[0]*(T_0-T_ref) + C_p[1]/2.0*(T_0-T_ref)**2 + C_p[2]/3.0*(T_0-T_ref)**3))
    c_1 = C_p[0]
    c_2 = C_p[1]/2.
    c_3 = C_p[2]/3.
    z = np.polynomial.polynomial.polyroots([c_0,c_1,c_2,c_3])
    z_r = z[z.imag==0][0].real
    T_1_analytic = T_ref + z_r

    assert 0.95 < (T_1_analytic-T_0)/(T_1_sim-T_0) < 1.05


def test_thermal_conductivity(res_path,tmp_path,np_rng,assert_allclose):
    def T_analytic(x, t, L, alpha, T_l, T_h):

        N = 10000
        sum_fourier = 0

        for i in range(1,N,1):

            omega = np.pi/L*(2*i-1)
            sum_fourier += (-1)**i/(2*i-1)*np.cos(omega*x)*np.exp(-alpha*t*omega**2)

        T = 0.5*(T_h+T_l) - 2./np.pi*(T_h-T_l)*sum_fourier
        return T

    grid = 'simple'
    load = 'no_stress'
    material = 'material'
    job = f'{grid}_{load}'

    K_11 = np_rng.random()+0.2
    C_p = np_rng.random()+0.2
    rho = np_rng.random()+0.2
    alpha = K_11/(rho*C_p)
    mat = damask.ConfigMaterial.load(res_path/f'{material}.yaml')
    del mat['phase']['matrix']['thermal']
    del mat['phase']['heatsource']['thermal']['K_33']
    del mat['phase']['heatsource']['thermal']['source']
    mat['phase']['heatsource']['lattice'] = 'cI'
    mat['phase']['heatsource']['rho'] = rho
    mat['phase']['heatsource']['thermal']['C_p'] = C_p
    mat['phase']['heatsource']['thermal']['K_11'] = K_11
    mat.save(tmp_path/f'{material}.yaml')

    unit_len = (np_rng.random()+1e-2)*1.0e-1
    T_ref = np_rng.integers(300,800)
    T_l = T_ref - np_rng.integers(50,200)
    T_h = T_ref + np_rng.integers(50,200)
    g = damask.GeomGrid(np.ones((32,2,2)),np.array([32,2,2])*unit_len)
    T_init = np.full((32,2,2), T_l)
    T_init[8:24,:,:] = T_h
    g.initial_conditions['T'] = T_init
    g.save(tmp_path/grid)

    l = damask.LoadcaseGrid.load(res_path/f'{load}.yaml')
    del l['loadstep'][0]['boundary_conditions']['mechanical']['P']
    l['loadstep'][0]['boundary_conditions']['mechanical']['F'] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    t = (g.size[0]**2/alpha)/8.
    l['loadstep'][0]['discretization']['t'] = t
    l['loadstep'][0]['discretization']['N'] = 200
    l.save(tmp_path/f'{load}.yaml')

    damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',wd=tmp_path)

    r = damask.Result(tmp_path/f'{job}.hdf5')
    T = [np.average(_.reshape((-1,2,2),order='F'),axis=(1,2)) for _ in r.get('T').values()]
    n = np_rng.integers(0, 200)
    T_1_sim = T[n][16:]
    x = np.linspace(0.5*unit_len,15.5*unit_len,16)
    T_1_analytic = np.array([T_analytic(x, n*t/200, 16*unit_len, alpha, T_l, T_h) for x in x])
    assert_allclose(T_1_analytic, T_1_sim, rtol=5e-2,atol=1.)

def test_thermal_dissipation(res_path,tmp_path,np_rng,assert_allclose):

    grid = 'simple'
    load = 'tensionX'
    material = 'material'
    job  = f'{grid}_{load}'
    config = 'config'

    l = damask.LoadcaseGrid.load(res_path/f'{load}.yaml')
    del l['loadstep'][0]
    t = np_rng.integers(20,50)
    l['loadstep'][0]['discretization']['t'] = t
    l['loadstep'][0]['discretization']['N'] = t
    l['loadstep'][0]['f_out'] = 1
    l.save(tmp_path/f'{load}.yaml')

    mat = damask.ConfigMaterial.load(res_path/f'{material}.yaml')
    del mat['phase']['matrix']
    del mat['material']
    del mat['phase']['heatsource']['thermal']['K_33']

    mat['phase']['heatsource']['mechanical']['output'] = ['F', 'P', 'S', 'L_p']
    mat['phase']['heatsource']['mechanical']['elastic'] = {'type':'Hooke', 'C_11': 243.9e+9, 'C_12': 157.0e+9, 'C_44': 117.9e+9}
    mat_config = damask.ConfigMaterial.load(res_path/f'{config}.yaml')
    mat['phase']['heatsource']['mechanical']['plastic'] = mat_config['phase']['IN625']['mechanical']['plastic']
    mat['phase']['heatsource']['thermal']['K_11'] = 0.0
    T_0 = np_rng.integers(300,350)
    C_p = (np_rng.random() + 0.2)*1.0e3
    rho = (np_rng.random() + 0.2)*1.0e4
    k = (np_rng.random() + 0.8)*5.0e-1

    mat['phase']['heatsource']['rho'] = rho
    mat['phase']['heatsource']['thermal']['C_p'] = C_p
    mat['phase']['heatsource']['thermal']['source'] = [{'type':'dissipation', 'kappa': k}]
    mat = mat.material_add(phase='heatsource',O=damask.Rotation.from_random(shape=8, rng_seed=np_rng), homogenization='direct')
    mat.save(tmp_path/f'{material}.yaml')

    g = damask.GeomGrid(np.ones((2,2,2)),np.ones(3)*1e-4)
    g.initial_conditions['T'] = T_0
    g.save(tmp_path/grid)

    damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',wd=tmp_path)
    r = damask.Result(tmp_path/f'{job}.hdf5')
    T_sim = np.array([np.average(_) for _ in r.get('T').values()])
    L_p = np.array([x.data for x in r.place('L_p').values()])
    S = np.array([x.data for x in r.place('S').values()])
    f_t = [np.average(np.sum(np.abs(_),axis=(1,2))) for _ in S*L_p]
    f_t_cum = np.cumsum(f_t)

    T_analytic = T_0 + f_t_cum*k/rho/C_p
    assert_allclose(T_analytic, T_sim, atol=0.5)
