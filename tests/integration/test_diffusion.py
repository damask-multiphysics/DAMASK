# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'diffusion'

def test_chemical_diffusion(res_path,tmp_path,np_rng,assert_allclose):
    def comp_analytic(x, t, L, D, C_l, C_h):

        N = 10000
        sum_fourier = 0

        for i in range(1,N,1):

            omega = np.pi/L*(2*i-1)
            sum_fourier += (-1)**i/(2*i-1)*np.cos(omega*x)*np.exp(-D*t*omega**2)

        comp = 0.5*(C_h + C_l) - 2./np.pi*(C_h - C_l)*sum_fourier
        return comp

    grid = 'diffusion_couple'
    load = 'no_deformation'
    material = 'material'
    job = f'{grid}_{load}'

    C_h = 0.2 + (np_rng.random()+0.5)*1.0e-1    # 0.25 < C_h < 0.35
    C_l = 0.2 - (np_rng.random()+0.5)*1.0e-1    # 0.05 < C_l < 0.15

    g = damask.GeomGrid.load(res_path/f'{grid}.vti')
    a = np.full((32,2,2), C_l)
    a[8:24,:,:] = C_h
    g.initial_conditions['Li'] = a
    g.initial_conditions['Al'] = 1.0 - a
    g.save(tmp_path/grid)

    mat = damask.ConfigMaterial.load(res_path/f'{material}.yaml')
    mat['phase']['Aluminum_1']['chemical']['components']['Li']['c_0'] = C_h
    mat['phase']['Aluminum_1']['chemical']['components']['Al']['c_0'] = 1 - C_h
    mat['phase']['Aluminum_2']['chemical']['components']['Li']['c_0'] = C_l
    mat['phase']['Aluminum_2']['chemical']['components']['Al']['c_0'] = 1 - C_l
    mat.save(tmp_path/f'{material}.yaml')

    D = 1.0
    l = damask.LoadcaseGrid.load(res_path/f'{load}.yaml')
    t = (g.size[0]**2/D)/8.
    l['loadstep'][0]['discretization']['t'] = t
    l['loadstep'][0]['discretization']['N'] = 500
    l.save(tmp_path/f'{load}.yaml')

    damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',wd=tmp_path)

    r = damask.Result(tmp_path/f'{job}.hdf5')
    Li_homog = np.array([x['homogenization'] for x in r.get('Li').values()])
    Li = [np.average(_.reshape((-1,2,2),order='F'),axis=(1,2)) for _ in Li_homog]

    tstep = np_rng.integers(50,400)
    time =  tstep*t/500
    x = np.linspace(0.5, 15.5, 16)
    C_analytic = np.array([comp_analytic(x, time, 16, D, C_l, C_h) for x in x])
    assert_allclose(C_analytic, Li[tstep][16:], rtol=5e-2,atol=1.)
