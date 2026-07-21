# SPDX-License-Identifier: AGPL-3.0-or-later
import os
import shutil

import pytest
import h5py
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'diffusion'

def h5py_dataset_iterator(g, prefix=''):
    for key,item in g.items():
        path = os.path.join(prefix, key)
        if isinstance(item, h5py.Dataset): # test for dataset
            yield (path, item)
        elif isinstance(item, h5py.Group): # test for group (go down)
            yield from h5py_dataset_iterator(item, path)


def h5py_compare_files(fname_reference,fname_current,assert_allclose):
    with h5py.File(fname_reference,'r') as ref, \
         h5py.File(fname_current,'r') as cur:
        for (path,dset) in h5py_dataset_iterator(ref):
            if path.startswith('setup'): continue
            if len(dset.dtype) > 0:
                for e in dset.dtype.fields:
                    assert np.all(dset[e]==cur[path][e])
            elif dset.dtype == 'int32':
                assert dset.size - np.count_nonzero(cur[path][()]==ref[path][()]) <= dset.size*0.0001
            else:
                if ref[path].size < 1: continue
                rtol = 0.025
                atol = max(np.abs(dset).max()*1e-5,1e-6)
                assert_allclose(cur[path],ref[path],rtol=rtol,atol=atol)


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

    damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',wd=tmp_path)

    r = damask.Result(tmp_path/f'{job}.hdf5')
    comp_homog = np.array([x for x in r.get('x').values()])
    Li_homog = np.array([e[:,0] for e in comp_homog])
    Li = [np.average(_.reshape((-1,2,2),order='F'),axis=(1,2)) for _ in Li_homog]

    tstep = np_rng.integers(50,400)
    time =  tstep*t/500
    x = np.linspace(0.5, 15.5, 16)
    C_analytic = np.array([comp_analytic(x, time, 16, D, C_l, C_h) for x in x])
    assert_allclose(C_analytic, Li[tstep][16:], rtol=5e-2,atol=1e-4)


def test_chemical_regularsolution_reference(res_path,tmp_path,assert_allclose,update):

    grid = 'regular'
    load = 'no_deformation'
    material = 'material_regular'
    job = f'{grid}_{load}'

    g = damask.GeomGrid.load(res_path/f'{grid}.vti')
    g.save(tmp_path/grid)

    mat = damask.ConfigMaterial.load(res_path/f'{material}.yaml')
    mat['phase']['phase_A']['chemical']['V_m'] = 0.0005212393458879053
    mat['phase']['phase_B']['chemical']['V_m'] = 0.0002967229062448962
    mat.save(tmp_path/f'{material}.yaml')

    l = damask.LoadcaseGrid.load(res_path/f'{load}.yaml')
    l['solver']['thermal'] = 'spectral'
    l['loadstep'][0]['discretization']['t'] = 1.5
    l['loadstep'][0]['discretization']['N'] = 100
    l.save(tmp_path/f'{load}.yaml')

    damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',wd=tmp_path)

    if update:
        shutil.copyfile(tmp_path/f'{job}.hdf5', res_path/f'regular_load_material.hdf5')

    h5py_compare_files(res_path/f'regular_load_material.hdf5', tmp_path/f'{job}.hdf5', assert_allclose)


def test_chemical_regularsolution_plausibility(res_path,tmp_path,np_rng,assert_allclose):

    grid = 'regular'
    load = 'no_deformation'
    material = 'material_regular'
    job = f'{grid}_{load}'

    g = damask.GeomGrid.load(res_path/f'{grid}.vti')
    g.save(tmp_path/grid)

    mat = damask.ConfigMaterial.load(res_path/f'{material}.yaml')
    V_m = (np_rng.random() + 10)*5.0e-5
    mat['phase']['phase_A']['chemical']['V_m'] = V_m
    V_m = (np_rng.random() + 5)*5.0e-5
    mat['phase']['phase_B']['chemical']['V_m'] = V_m
    mat.save(tmp_path/f'{material}.yaml')

    t = np_rng.uniform(1.2,1.6)
    N = np_rng.integers(95,105)
    l = damask.LoadcaseGrid.load(res_path/f'{load}.yaml')
    l['solver']['thermal'] = 'spectral'
    l['loadstep'][0]['discretization']['t'] = t
    l['loadstep'][0]['discretization']['N'] = N
    l.save(tmp_path/f'{load}.yaml')

    damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',wd=tmp_path)

    # check 1: decrease in chemical potential gradient (equivalently difference) over time
    r = damask.Result(tmp_path/f'{job}.hdf5')
    mu_Nb = np.array([np.average(mu[:,0].reshape(g.cells,order='F'),axis=(1,2)) for mu in r.get('mu').values()])
    delta_mu_Nb = mu_Nb.max(axis=1) - mu_Nb.min(axis=1)
    assert (delta_mu_Nb[:-1]-delta_mu_Nb[1:] > 0).all()

    # check 2: symmetric result for a symmetric initial profile
    inc = np_rng.integers(1,N)
    x_Nb = np.average(r.view(increments=inc).get('x')[:,0].reshape(g.cells,order='F'),axis=(1,2))
    assert_allclose(x_Nb[0:g.cells[0]//2], x_Nb[g.cells[0]:g.cells[0]//2-1:-1], rtol=5e-2,atol=1e-4)
