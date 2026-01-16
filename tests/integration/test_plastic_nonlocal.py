# SPDX-License-Identifier: AGPL-3.0-or-later
import itertools

import pytest
import numpy as np
from scipy import stats

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'plastic_nonlocal'


@pytest.mark.parametrize('integrator',['FPI','Euler','Euler_Adaptive','RK4','RKCK45'])
def test_density_conservation(res_path,tmp_path,copy_files,integrator):

   load='shearXY'
   grid='rim'
   material = 'material'
   job=f'{grid}_{load}'
   nominal_rho = 2.25e14

   copy_files(res_path,tmp_path,[f'{load}.yaml',f'{grid}.vti',f'{material}.yaml'])

   config_numerics = damask.YAML.load(res_path/'numerics.yaml')
   config_numerics['phase']['integrator_state'] = integrator
   config_numerics.save(tmp_path/'numerics.yaml')

   stdout,stderr = damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml '+
                                   f'-n numerics.yaml -j {job}',wd=tmp_path)
   print(f'{stdout}\n{stderr}')

   result = damask.Result(tmp_path/f'{job}.hdf5')
   result.add_calculation('np.abs(#rho_u_ed_pos#)+np.abs(#rho_u_ed_neg#)+'+\
                          'np.abs(#rho_u_sc_pos#)+np.abs(#rho_u_sc_neg#)+'+\
                          'np.abs(#rho_b_ed_pos#)+np.abs(#rho_b_ed_neg#)+'+\
                          'np.abs(#rho_b_sc_pos#)+np.abs(#rho_b_sc_neg#)+'+\
                          '#rho_d_ed#+#rho_d_sc#',
                          'rho','1/m²','total dislocation density')
   rho = np.array([np.average(rho_field,axis=0) for rho_field in result.place('rho').values()])

   print(f'{np.min(rho)} < rho < {np.max(rho)}')
   assert (nominal_rho*.999 <rho).all() and (rho < nominal_rho*1.001).all()

@pytest.mark.parametrize('flux',[True,False])
def test_Hall_Petch(res_path,tmp_path,copy_files,assert_allclose,np_rng,flux):
   load = 'tensionY'
   copy_files(res_path,tmp_path,[f'{load}.yaml'])

   N_grains = 18
   cells = np.array([16]*3)
   seeds = damask.seeds.from_random(size=np.ones(3),cells=cells,N_seeds=N_grains,rng_seed=np_rng)
   grid = damask.GeomGrid.from_Voronoi_tessellation(size=np.ones(3),cells=cells,seeds=seeds)

   config_material = damask.ConfigMaterial.load(res_path/'material.yaml')
   del(config_material['phase']['Aluminum_empty'])
   del(config_material['material'])
   config_material['phase']['Aluminum']['mechanical']['plastic']['flux'] = flux

   config_material.material_add(phase='Aluminum',homogenization='SX',
                                O = damask.Rotation.from_random(N_grains,rng_seed=np_rng)) \
                               .save(tmp_path/'material.yaml')

   stress_prev = None
   for size in [1e-2,1e-3,5e-4,1e-4]:
      grid.size = np.array([size]*3)
      grid.save(tmp_path/f'{size}.vti')

      damask.util.run(f'DAMASK_grid -l {load}.yaml -g {size}.vti -m material.yaml --wd {tmp_path} -j {size}_{load}')
      r = damask.Result(tmp_path/f'{size}_{load}.hdf5')
      stress = np.array([np.average(P[:,1,1]) for P in r.view(increments=[40,50,60]).get('P').values()])

      if flux:     assert stress_prev is None or np.all(stress > stress_prev)
      if not flux and stress_prev is not None:
         assert_allclose(stress,stress_prev)

      stress_prev = stress

def test_density_initialization(res_path,tmp_path,np_rng):

   load='none'
   grid='simple'
   material = 'material'

   rho_u = {}
   for c,s in itertools.product(['ed','sc'],['pos','neg']):
      rho_u[f'{c}_{s}'] = 10**np_rng.uniform(13,15)
   sigma_rho_u = 10**np_rng.uniform(13,15)

   load_case = damask.YAML(solver={'mechanical':'spectral_basic'},loadstep=[])
   load_case['loadstep'].append({'boundary_conditions':{'mechanical':{'F':np.eye(3).tolist()}},
                                 'discretization':{'t':1.,'N':1}})
   load_case.save(tmp_path/f'{load}.yaml')

   g = damask.GeomGrid(np.zeros(32**3).reshape([32]*3).astype(int),np.ones(3)*1e-3)
   g.save(tmp_path/grid)

   config_material = damask.ConfigMaterial.load(res_path/f'{material}.yaml')
   config_material['phase']['Aluminum']['mechanical']['plastic']['N_slip'] = [12]
   config_material['phase']['Aluminum']['mechanical']['plastic']['sigma_rho_u'] = sigma_rho_u
   for k,v in rho_u.items():
      config_material['phase']['Aluminum']['mechanical']['plastic'][f'rho_u_{k}_0'] = [v]
   config_material.save(tmp_path/f'{material}.yaml')


   stdout,stderr = damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {grid}_{load}',
                                   wd=tmp_path)

   r = damask.Result(tmp_path/f'{grid}_{load}.hdf5').view(increments=0)
   passes = []
   for k,v in rho_u.items():
      rho = r.get(f'rho_u_{k}')
      passes.append(abs(np.std(rho)-sigma_rho_u)/sigma_rho_u < 0.1)
      passes.append(abs(np.mean(rho)-v)/v < 0.1)
      passes.append(stats.normaltest(rho,axis=None)[1]>0.01)

   assert np.count_nonzero(passes) > 8
