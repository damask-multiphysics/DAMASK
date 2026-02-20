# SPDX-License-Identifier: AGPL-3.0-or-later
from pathlib import Path
import itertools

import pytest
import numpy as np
import h5py

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'configexamples'


config_homogenization = \
   [
    # mechanical
    ('direct','pass_direct', None,None),
    ('bicrystal','isostrain_polycrystal', None,None),
    ('8grains','RGC_8grains', None,None),
    # mechanical + thermal
    ('direct','pass_direct','pass_direct',None),
    ('8grains','isostrain_polycrystal','isotemperature_polycrystal',None),
    # mechanical + thermal
    ('direct','pass_direct','pass_direct','pass_direct'),
   ]


config_phase = \
   [
    # elastic
    ('fcc', 'Hooke_TWIP-steel', None, None, None, None, None, None),
    ('bcc', 'Hooke_Nb', None, None, None, None, None, None),
    ('bcc', 'Hooke_Ta', None, None, None, None, None, None),
    ('bct', 'Hooke_Sn-beta', None, None, None, None, None, None),
    ('hcp', 'Hooke_SiC-6H', None, None, None, None, None, None),
    ('Al', 'Hooke_Al', None, None, None, None, None, None),
    ('Cd', 'Hooke_Cd', None, None, None, None, None, None),
    ('Fe', 'Hooke_C50E-martensite', None, None, None, None, None, None),
    ('Mg', 'Hooke_Mg', None, None, None, None, None, None),
    ('Ni', 'Hooke_Ni', None, None, None, None, None, None),
    ('Ni', 'Hooke_IN625', None, None, None, None, None, None),
    ('Sn-beta', 'Hooke_Sn-beta', None, None, None, None, None, None),
    ('Pt', 'Hooke_Pt', None, None, None, None, None, None),
    ('Si', 'Hooke_Si', None, None, None, None, None, None),
    ('Ti', 'Hooke_Ti', None, None, None, None, None, None),
    ('W', 'Hooke_W', None, None, None, None, None, None),
    ('W', 'Hooke_W-3at.%Re', None, None, None, None, None, None),
    ('W', 'Hooke_W-10at.%Re', None, None, None, None, None, None),
    # elastic + damage
    ('Al', 'Hooke_Al', None, None, 'anisobrittle_cubic', None, None, None),
    ('Fe', 'Hooke_Fe', None, None, 'isobrittle_generic', None, None, None),
    # elastic + plastic
    ('fcc', 'Hooke_vanishing-Poisson-ratio', 'dislotwin_alpha-Brass-shearbanding', None, None, None, None, None),
    ('fcc', 'Hooke_vanishing-Poisson-ratio', 'isotropic_free-surface', None, None, None, None, None),
    ('fcc', 'Hooke_TWIP-steel', 'dislotwin_TWIP-TRIP', None, None, None, None, None),
    ('Al', 'Hooke_Al', 'phenopowerlaw_Al', None, None, None, None, None),
    ('Al', 'Hooke_Al', 'phenopowerlaw_AA6022-T4', None, None, None, None, None),
    ('Al', 'Hooke_Al', 'nonlocal_Al', None, None, None, None, None),
    ('Au', 'Hooke_Au', 'phenopowerlaw_Au', None, None, None, None, None),
    ('Cu', 'Hooke_Cu', 'phenopowerlaw_Cu', None, None, None, None, None),
    ('Fe', 'Hooke_Fe', 'phenopowerlaw_bcc-martensite', None, None, None, None, None),
    ('Fe', 'Hooke_Fe', 'phenopowerlaw_DP-steel-ferrite', None, None, None, None, None),
    ('Fe', 'Hooke_Fe', 'phenopowerlaw_polygonal-ferrite', None, None, None, None, None),
    ('Fe', 'Hooke_Fe', 'dislotwin_IF-steel', None, None, None, None, None),
    ('Mg', 'Hooke_Mg', 'phenopowerlaw_Mg', None, None, None, None, None),
    ('Ni', 'Hooke_Ni', 'nonlocal_Ni', None, None, None, None, None),
    ('Ni', 'Hooke_IN625', 'phenopowerlaw_IN625', None, None, None, None, None),
    ('Ni', 'Hooke_IN625', 'dislotwin_IN625', None, None, None, None, None),
    ('Pt', 'Hooke_Pt', 'phenopowerlaw_Pt-5%Cu', None, None, None, None, None),
    ('Sn-beta', 'Hooke_Sn-beta', 'phenopowerlaw_Sn-beta', None, None, None, None, None),
    ('Ti', 'Hooke_Ti', 'phenopowerlaw_Ti', None, None, None, None, None),
    ('W', 'Hooke_W', 'dislotungsten_W', None, None, None, None, None),
    ('X5CrNi18-10', 'Hooke_X5CrNi18-10', 'kinehardening_X2CrNiMo18-15-4', None, None, None, None, None),
    # elastic + plastic + damage
    ('Cu', 'Hooke_Cu', 'phenopowerlaw_Cu', None, 'isobrittle_generic', None, None, None),
    ('Fe', 'Hooke_Fe', 'dislotwin_IF-steel', None, 'anisobrittle_cubic', None, None, None),
    # elastic + thermal
    ('Al', 'Hooke_Al', None, None, None, 'Al', None, None),
    ('Ag', 'Hooke_Ag', None, None, None, 'Ag', None, None),
    ('Pt', 'Hooke_Pt', None, None, None, 'Pt', None, None),
    # elastic + thermal + source
    ('Fe', 'Hooke_Fe', None, None, None, 'Fe', 'externalheat_ramp-and-hold', None),
    ('Ni', 'Hooke_Ni', None, None, None, 'Ni', 'externalheat_ramp-and-hold', None),
    ('Ti', 'Hooke_Ti', None, None, None, 'adiabatic', 'externalheat_ramp-and-hold', None),
    # elastic + plastic + thermal + source
    ('Al', 'Hooke_Al', 'phenopowerlaw_Al', None, None, 'Al', 'dissipation_generic', None),
    ('Au', 'Hooke_Au', 'phenopowerlaw_Au', None, None, 'fast-convection', 'externalheat_ramp-and-hold', None),
    ('Au', 'Hooke_Au', 'phenopowerlaw_Au', None, None, 'adiabatic', 'dissipation_generic', None),
    ('Fe', 'Hooke_Fe', 'dislotwin_IF-steel', None, None, 'steel-0.5C', 'dissipation_generic', None),
    # elastic + eigen
    ('bcc', 'Hooke_W', None, 'thermalexpansion_SiC-beta', None, None, None, None),
    ('hcp', 'Hooke_Mg', None, 'thermalexpansion_SiC-alpha', None, None, None, None),
    ('Mg', 'Hooke_Mg', None, 'thermalexpansion_Mg', None, None, None, None),
    ('W', 'Hooke_W', None, 'thermalexpansion_W', None, None, None, None),
    # elastic + eigen + thermal + source
    ('Al', 'Hooke_Al', None, 'thermalexpansion_Al', None, 'Al', 'externalheat_ramp-and-hold', None),
    ('Au', 'Hooke_Au', None, 'thermalexpansion_Au', None, 'Au', 'externalheat_ramp-and-hold', None),
    ('Cu', 'Hooke_Cu', None, 'thermalexpansion_Cu', None, 'Cu', 'externalheat_ramp-and-hold', None),
    ('Fe', 'Hooke_Fe', None, 'thermalexpansion_Fe', None, 'steel-0.5C', 'externalheat_ramp-and-hold', None),
    ('Fe', 'Hooke_Fe', None, 'thermalexpansion_C35E', None, 'steel-0.5C', 'externalheat_ramp-and-hold', None),
    ('Fe', 'Hooke_Fe', None, 'thermalexpansion_X20Cr13', None, 'steel-0.5C', 'externalheat_ramp-and-hold', None),
    ('Si', 'Hooke_Si', None, 'thermalexpansion_Si', None, 'Si','externalheat_ramp-and-hold', None),
    ('Sn-beta', 'Hooke_Sn-beta', None, 'thermalexpansion_Sn-beta', None, 'Sn-beta','externalheat_ramp-and-hold', None),
    ('W', 'Hooke_W', None, 'thermalexpansion_W', None, 'W', 'externalheat_ramp-and-hold', None),
    ('X2CrNi18-10', 'Hooke_X6CrNiMo17-12-2', None, 'thermalexpansion_X2CrNi18-10', None, 'X2CrNi18-10',
                                                                                     'externalheat_ramp-and-hold', None),
    ('X2CrNiMo17-12-2', 'Hooke_X2CrNiMo17-12-2', None, 'thermalexpansion_X2CrNiMo17-12-2', None, 'X2CrNiMo17-12-2',
                                                                                     'externalheat_ramp-and-hold', None),
    ('X5CrNi18-10', 'Hooke_X5CrNi18-10', None, 'thermalexpansion_X5CrNi18-10', None, 'X5CrNi18-10',
                                                                                     'externalheat_ramp-and-hold', None),
    # elastic + chemical
    ('Al', 'Hooke_Al', None, None, None, None, None, 'quadEnergy_example'),
   ]


def test_coverage_homogenization(damask_root):
    existing = set((damask_root/'examples/config/homogenization').glob('**/*.yaml'))

    tested   = set()
    path = ['.','mechanical','thermal','damage']
    for combination in config_homogenization:
        tested.update([damask_root/f'examples/config/homogenization/{path[i]}/{config}.yaml'
                       for i,config in enumerate(combination) if config is not None])

    assert len(existing ^ tested) == 0, existing ^ tested


@pytest.mark.parametrize('N_constituents,mechanical,thermal,damage',config_homogenization)
def test_homogenization(damask_root,tmp_path,
                        N_constituents,mechanical,thermal,damage):
    def get_output(config,label='DAMASK'):
        output = {}
        homogenization = config['homogenization'][label]
        for t in ['mechanical','thermal','damage']:
            if t in homogenization.keys():
                output[t] = homogenization[t].get('output', [])
        return output

    grid = 'SX'
    load = 'none'
    material = 'material'
    label = 'DAMASK'
    base_path = Path(damask_root/'examples/config/homogenization')

    g = damask.GeomGrid(np.zeros([2,1,1]),np.ones(3))

    load_case = damask.YAML(solver={'mechanical':'spectral_basic'},loadstep=[])
    load_case['loadstep'].append({'boundary_conditions':{'mechanical':{'F':np.eye(3).tolist()}},
                                  'discretization':{'t':1.,'N':1}})

    config = damask.ConfigMaterial()
    config['phase']['iso'] = {'mechanical': {'elastic':{'type':'Hooke','C_11':100,'C_12':50,'C_44':25}},
                              'lattice': 'cF',
                              'rho': 1}

    config['homogenization'][label] = damask.YAML.load(base_path/f'{N_constituents}.yaml')
    N = int(config['homogenization'][label]['N_constituents'])

    config['homogenization'][label]['mechanical'] \
        = damask.YAML.load(base_path/'mechanical'/f'{mechanical}.yaml')

    if thermal is not None:
        config['homogenization'][label]['thermal'] \
            = damask.YAML.load(base_path/'thermal'/f'{thermal}.yaml')
        config['phase']['iso']['thermal'] = {'C_p':1, 'K_11':1}
        g.initial_conditions['T'] = 200.0
        load_case['solver']['thermal'] = 'spectral'

    if damage is not None:
        config['homogenization'][label]['damage'] \
            = damask.YAML.load(base_path/'damage'/f'{damage}.yaml')
        config['phase']['iso']['damage'] \
            = damask.YAML.load(damask_root/'examples'/'config'/'phase'/'damage'/'isobrittle_generic.yaml')
        load_case['solver']['damage'] = 'spectral'
        g.initial_conditions['phi'] = 1.0



    config.material_add(phase='iso',O=damask.Rotation().broadcast_to((1,N)),                    # noqa
                        homogenization=label).save(tmp_path/f'{material}.yaml')
    load_case.save(tmp_path/f'{load}.yaml')
    g.save(tmp_path/f'{grid}.vti')

    out,err = damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {grid}_{load}',
                              wd=tmp_path)
    print(out)
    print(err)

    output_expected = get_output(config)
    r = damask.Result(tmp_path/f'{grid}_{load}.hdf5').view(increments=0,phases=False)
    output = r.get(list(itertools.chain(*output_expected.values())),flatten=False,prune=False)
    for k,v in output_expected.items():
        assert set(output['increment_0']['homogenization'][label][k].keys()) == set(v)

def test_coverage_phase(damask_root):
    existing = set((damask_root/'examples/config/phase').glob('**/*.yaml'))

    tested   = set()
    path = ['.','mechanical/elastic','mechanical/plastic','mechanical/eigen',
            'damage','thermal','thermal/source', 'chemical',
           ]
    for combination in config_phase:
        tested.update([damask_root/f'examples/config/phase/{path[i]}/{config}.yaml'
                       for i,config in enumerate(combination) if config is not None])

    assert len(existing ^ tested) == 0, existing ^ tested


@pytest.mark.parametrize('lattice,elastic,plastic,eigen,damage,thermal,source,chemical',config_phase)
def test_phase(damask_root,tmp_path,
               lattice,elastic,plastic,eigen,damage,thermal,source,chemical):

    def get_output(config,label='DAMASK'):
        output = {}
        phase = config['phase'][label]
        for t in ['mechanical','thermal','damage','chemical']:
            if t in phase.keys():
                output[t] = []
                if t == 'mechanical':
                    for m in ['elastic','plastic']:
                        if m in phase[t].keys():
                            output[t] += phase[t][m].get('output', [])
                    if 'eigen' in phase[t].keys():
                        for e in phase[t]['eigen']:
                            output[t] += e.get('output', [])
                else:
                    output[t] += phase[t].get('output', [])
        return output

    grid = 'SX'
    load = 'none'
    material='material'
    label = 'DAMASK'
    base_path = Path(damask_root/'examples/config/phase')

    g = damask.GeomGrid(np.zeros([2,1,1]),np.ones(3))

    load_case = damask.YAML(solver={'mechanical':'spectral_basic'},loadstep=[])
    load_case['loadstep'].append({'boundary_conditions':{'mechanical':{'F':np.eye(3).tolist()}},
                                  'discretization':{'t':1.,'N':1}})

    config = damask.ConfigMaterial()
    config['homogenization']['SX'] = {'N_constituents':1,'mechanical':{'type':'pass'}}

    config['phase'][label] = damask.YAML.load(base_path/f'{lattice}.yaml')
    config['phase'][label]['mechanical'] = {}
    config['phase'][label]['mechanical']['elastic'] \
        = damask.YAML.load(base_path/'mechanical'/'elastic'/f'{elastic}.yaml')

    if plastic is not None:
        config['phase'][label]['mechanical']['plastic'] \
            = damask.YAML.load(base_path/'mechanical'/'plastic'/f'{plastic}.yaml')

    if eigen is not None:
        config['phase'][label]['mechanical']['eigen'] \
            = [damask.YAML.load(base_path/'mechanical'/'eigen'/f'{eigen}.yaml')]

    if thermal is not None:
        config['phase'][label]['thermal'] \
            = damask.YAML.load(base_path/'thermal'/f'{thermal}.yaml')
        config['homogenization']['SX']['thermal'] = {'type':'pass'}
        g.initial_conditions['T'] = 200.0
        load_case['solver']['thermal'] = 'spectral'

    if source is not None:
        config['phase'][label]['thermal']['source'] \
            = [damask.YAML.load(base_path/'thermal'/'source'/f'{source}.yaml')]

    if damage is not None:
        config['phase'][label]['damage'] \
            = damask.YAML.load(base_path/'damage'/f'{damage}.yaml')
        config['homogenization']['SX']['damage'] = {'type':'pass'}
        load_case['solver']['damage'] = 'spectral'
        g.initial_conditions['phi'] = 1.0

    config.material_add(phase=label,O=[1.,0.,0.,0.],homogenization='SX').save(tmp_path/f'{material}.yaml') # noqa
    load_case.save(tmp_path/f'{load}.yaml')
    g.save(tmp_path/f'{grid}.vti')

    out, err = damask.util.run(f'damask_grid -l {load}.yaml -g {grid}.vti -m {material}.yaml -j {grid}_{load}',
                               wd=tmp_path)
    print(out)
    print(err)

    output_expected = get_output(config)
    r = damask.Result(tmp_path/f'{grid}_{load}.hdf5').view(increments=0,homogenizations=False)
    output = r.get(list(itertools.chain(*output_expected.values())),flatten=False,prune=False)
    for k,v in output_expected.items():
        assert set(output['increment_0']['phase'][label][k].keys()) == set(v)

    with open(tmp_path/f'{material}.yaml') as f_ASCII, h5py.File(tmp_path/f'{grid}_{load}.hdf5','r') as f_HDF5:
        assert f_HDF5['setup/material.yaml'][0].decode() == f_ASCII.read()
