#!/usr/bin/env python3

from pathlib import Path
import os

import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt

import damask


# material.yaml
example_dir = Path.home()/'DAMASK/examples/config/phase'

phase = damask.ConfigMaterial.load(example_dir/'Cu.yaml')
phase['mechanical'] = {'output': ['F','P']}
phase['mechanical']['elastic'] = damask.ConfigMaterial.load(example_dir/'mechanical/elastic/Hooke_Cu.yaml')
phase['mechanical']['plastic'] = damask.ConfigMaterial.load(example_dir/'mechanical/plastic/phenopowerlaw_Cu.yaml')

mat = damask.ConfigMaterial()
mat['phase']['Cu'] = phase
mat['homogenization']['SX'] = {'N_constituents': 1,'mechanical': {'type': 'pass'}}

# load
load = damask.Config()
load['solver'] = {'mechanical': 'spectral_basic'}
load['loadstep'] = [{'boundary_conditions': {'mechanical': {'dot_F': [[3e-3,'x','x'],[0,'x','x'],[0,0,'x']],
                                                            'P': [['x',0,0],['x',0,0],['x','x',0]]}},
                           'discretization': {'t': 100, 'N': 100}}]

# grid
grid = damask.Grid(np.zeros([2,2,2],int),np.ones(3)*1e-5)


# crystallographic directions
samples = {'001':[[0,0,1],[0,1,0]],'111':[[1,1,1],[0,1,-1]]}

# reference results
reference = {}
for label,directions in samples.items():
    ref_data = damask.Table.load('Takeuchi1975_Fig3b_'+label+'.txt')
    reference[label] = [ref_data.get('epsilon')[:,0]*.01,ref_data.get('sigma')[:,0]*1e6]

results = {}
for label,directions in samples.items():
    os.makedirs(label, exist_ok=True)
    
    grid.save('/'.join([label,'2x2x2.vti']))
    load.save('/'.join([label,'tensionX.yaml']))
    
    ori = damask.Orientation.from_directions(directions[0],directions[1],lattice='cF')
    mat_c = mat.material_add(O=ori,phase='Cu',homogenization='SX')
    mat_c.save('/'.join([label,'material.yaml']))
    
    damask.util.execute(f'DAMASK_grid -l tensionX.yaml -g 2x2x2.vti --wd {label}')
    
    r = damask.Result('/'.join([label,'2x2x2_tensionX.hdf5']))
    results[label] = [np.array([np.average(F[:,0,0]-1) for F in r.get('F').values()]),
                      np.array([np.average(P[:,0,0]) for P in r.get('P').values()])]

    color = 'tab:blue' if label == '111' else 'tab:orange'
    plt.plot(results[label][0],results[label][1],label=label+'_sim',color=color,linestyle=':')
    plt.plot(reference[label][0],reference[label][1],label=label+'_exp',color=color,linestyle='-')
    plt.savefig('Takeuchi1975_Fig3b.pdf')
