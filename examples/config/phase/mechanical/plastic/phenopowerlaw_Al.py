#!/usr/bin/env python3

from pathlib import Path
import os

import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt

import damask


# material.yaml
example_dir = Path.home()/'DAMASK/examples/config/phase'

phase = damask.ConfigMaterial.load(example_dir/'Al.yaml')
phase['mechanical'] = {'output': ['F','P']}
phase['mechanical']['elastic'] = damask.ConfigMaterial.load(example_dir/'mechanical/elastic/Hooke_Al.yaml')
phase['mechanical']['plastic'] = damask.ConfigMaterial.load(example_dir/'mechanical/plastic/phenopowerlaw_Al.yaml')

mat = damask.ConfigMaterial()
mat['phase']['Al'] = phase
mat['homogenization']['SX'] = {'N_constituents': 1,'mechanical': {'type': 'pass'}}

# load
load = damask.Config()
load['solver'] = {'mechanical': 'spectral_basic'}
load['loadstep'] = [{'boundary_conditions': {'mechanical': {'dot_F': [[4.5e-3,'x','x'],[0,'x','x'],[0,0,'x']],
                                                            'P': [['x',0,0],['x',0,0],['x','x',0]]}},
                           'discretization': {'t': 65, 'N': 100}}]

# grid
grid = damask.Grid(np.zeros([2,2,2],int),np.ones(3)*1e-5)

# crystallographic directions
samples = {'001':[[0,0,1],[0,1,0]],'111':[[1,1,1],[0,1,-1]]}

# reference results
reference = {}
for label,directions in samples.items():
    ref_data = damask.Table.load('HosfordEtAl1960_Fig5_'+label+'.txt')
    reference[label] = [ref_data.get('epsilon')[:,0],ref_data.get('sigma')[:,0]*6894.76e3] # kpsi to Pa

results = {}
for label,directions in samples.items():
    
    os.makedirs(label, exist_ok=True)

    grid.save('/'.join([label,'2x2x2.vti']))
    load.save('/'.join([label,'tensionX.yaml']))

    ori = damask.Orientation.from_directions(directions[0],directions[1],lattice='cF')
    mat_c = mat.material_add(O=ori,phase='Al',homogenization='SX')
    mat_c.save('/'.join([label,'material.yaml']))
    
    damask.util.execute(f'DAMASK_grid -l tensionX.yaml -g 2x2x2.vti --wd {label}')
    
    r = damask.Result('/'.join([label,'2x2x2_tensionX.hdf5']))
    r.add_strain()
    r.add_stress_Cauchy()
    results[label] = [np.array([np.average(epsilon[:,0,0]) for epsilon in r.get('epsilon_V^0.0(F)').values()]),
                      np.array([np.average(sigma[:,0,0]) for sigma in r.get('sigma').values()])]


    color = 'tab:blue' if label == '111' else 'tab:orange'
    plt.plot(results[label][0],results[label][1],label=label+'_sim',color=color,linestyle=':')
    plt.plot(reference[label][0],reference[label][1],label=label+'_exp',color=color,linestyle='-')
    plt.savefig('HosfordEtAl1960_Fig5.pdf')
