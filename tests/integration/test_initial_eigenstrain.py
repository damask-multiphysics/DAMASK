# SPDX-License-Identifier: AGPL-3.0-or-later
import numpy as np

import damask

def test_initial_eigenstrain(tmp_path,assert_allclose,np_rng):

      grid = 'simple'
      load = 'no_stress'
      material = 'material'
      job  = f'{grid}_{load}'

      g = np.array([2,3,4])
      N = np.prod(g)

      O_0 = damask.Rotation.from_random(N,rng_seed=np_rng)
      V_e_0 = damask.tensor.symmetric(np.broadcast_to(np.eye(3),(N,3,3))+ 0.05*(np_rng.random((N,3,3))-0.5))

      damask.GeomGrid(np.arange(N,dtype=int).reshape(g,order='F'),
                  np.ones(3)*1e-6)\
            .save(tmp_path/grid)

      damask.ConfigMaterial({
            'homogenization': {'SX':{'N_constituents':1,'mechanical':{'type':'pass'}}},
            'phase': {'A':{'lattice':'cF','mechanical':{'output':['O','F','F_i','F_e','F_p'],
                                                        'elastic': {'type':'Hooke',
                                                                    'C_11':105.0e+9,
                                                                    'C_12':65.0e+9,
                                                                    'C_44':30.0e+9}}}}
                        })\
            .material_add(phase='A',
                          homogenization='SX',
                          V_e=V_e_0,
                          O=O_0)\
            .save(tmp_path/f'{material}.yaml')

      F = [['x', 0 , 0 ],
           [ 0, 'x', 0 ],
           [ 0,  0 ,'x']]
      P = [[ 0, 'x','x'],
           ['x', 0, 'x'],
           ['x','x', 0 ]]

      damask.YAML(solver={'mechanical':'spectral_basic'},
                  loadstep=[{'boundary_conditions':{'mechanical':{'F':F,'P':P}},
                             'discretization':{'t':1.,'N':1}}])\
            .save(tmp_path/f'{load}.yaml')

      damask.util.run(f'DAMASK_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {job}',wd=tmp_path)
      r = damask.Result(tmp_path/f'{job}.hdf5').view(increments=0)

      F_e = r.get('F_e')
      F_i = r.get('F_i')
      F_p = r.get('F_p')
      O   = damask.Rotation(r.get('O'))

      assert_allclose(O@np.linalg.inv(V_e_0),F_i)
      assert_allclose(O_0,O)
      assert_allclose(np.eye(3),F_e@F_i@F_p)
      assert_allclose(V_e_0,damask.mechanics.stretch_left(F_e))
