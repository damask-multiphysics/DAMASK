import pytest
import numpy as np
import damask

""" ----------------------------------------------------------------------------

    *** TEST COHERENCE OF BOUNDARY CONDITIONS ***

    Check boundary conditions are properly applied, such that:
    - Discretization-independency i.e. the number of discretization steps `N`
      does not change the solution.
      - Test versions 1, 3, 4
    - Rates and aim (target) equivalence i.e. the solution is the same when BC
      are applied using equivalent rates (x_dot) or aim (x) values
      - Example: {x_dot=5, t = 2} is equivalent to {x = 10, t = []}
    - Time-independence for aim/target BC (version 2 below)
      - Example: {x = 10, t = ..} should give the same result regardless of 't'
      - Test versions 2
    - For the linear case, additive BC
      - Test versions 3, 4

    Parametrization:
     - Dimension (2 or 3)
     - Boundary condition type (u_dot)

    Randomization:
     - Elastic constants (C_11, C_12, ...)
     - Orientations (in the material.yaml file)
     - Magnitude of the load boundary condition
     - Load case configuration:
        - Number of steps
        - Number of discretizations (N) per step

    Resources (in resources/mesh_BC_coherence):
    - Load (-l)
      - load_{2/3}D.yaml
    - Material (-m)
      - material.yaml
    - Mesh (-g)
      - mesh_{2/3}D.msh (single QUAD/HEXA element)
    - Numerics configuration (-n)
      - numerics.yaml

    Tests
    - Versions 0..4:
        Checks final displacement for the same (total) load applied using a
        different number of load steps, discretization (N) and time step
        length (t):
        - Version 0: single step, N = 1, t = 1.0 (reference)
        - Version 1: single step, N > 1, t = 1.0
        - Version 2: single step, N = 1, t > 1.0
        - Version 3: k steps, N >= 1, t = 1.0
        - Version 4: 2 steps, N >= 1, t = 1.0 (first and last), plus k steps
                     in between with free BC ([x,x,x])
    - Version 5:
        Checks the final displacement is zero after applying a "cyclic" load,
        simulated by consecutive steps of (total) force {+F, -F, +F, -F...}.

    - All versions (0..5) are run for both rates (f_dot) and aim (f)
      boundary conditions. Version 0 with rate f_dot is used as reference.

    ------------------------------------------------------------------------ """


@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'mesh_BC_coherence'


def bc_setup(np_rng, n_D, BC_type):
    """Randomized displacement boundary conditions."""
    # Set BC (along Y in 2D, along Z in 3D)
    scale = 1.0e+6 if 'BC_type' == 'f' else 1.0e-2
    bc_val = ['x', 'x', 'x']
    bc_val[n_D-1] = scale*np_rng.uniform(0.5,1)

    return bc_val


def load_dict(n_D, bc_val, bc_key, N = 1, t = 1.0):
    """Create a new load dictionary."""
    if n_D == 2:
        mech = [{'label': 'load',   bc_key:  bc_val},
                {'label': 'fix_uv', 'u_dot': [0.0, 0.0, 'x']},
                {'label': 'fix_v',  'u_dot': ['x', 0.0, 'x']}]
    else:
        mech = [{'label': 'load',    bc_key:  bc_val},
                {'label': 'fix_uvw', 'u_dot': [0.0, 0.0, 'x']},
                {'label': 'fix_u',   'u_dot': [0.0, 'x', 'x']},
                {'label': 'fix_w',   'u_dot': ['x', 'x', 0.0]}]

    d = {'boundary_conditions': {'mechanical': mech},
         'discretization': {'t': t, 'N': N},
         'f_out': 1
        }

    return d


def load_setup(load_config, np_rng, bc_val, n_D, key, version):
    """Randomized displacement boundary conditions."""
    # Randomize the number of steps and discretizations per step
    if version == 0:
        load_config['loadstep'][0] = load_dict(n_D, bc_val, key)
    elif version == 1:
        load_config['loadstep'][0] = load_dict(n_D, bc_val, key, N = np_rng.integers(1,10))
    elif version == 2:
        k = np_rng.integers(2,10)
        bc_new = bc_val[:]
        if key[-4:] == '_dot': bc_new[n_D-1] /= k
        load_config['loadstep'][0] = load_dict(n_D, bc_new, key, N = np_rng.integers(1,10), t = k)
    elif version == 3:
        k = np_rng.integers(2,6)
        bc_new = bc_val[:]
        bc_new[n_D-1] /= k
        load_config['loadstep'][0] = load_dict(n_D, bc_new, key, N = np_rng.integers(1,10))
        for m in np.arange(1,k):
            if not key[-4:] == '_dot': bc_new[n_D-1] += bc_val[n_D-1]/k
            load_config['loadstep'].append(load_dict(n_D, bc_new, key, N = np_rng.integers(1,10)))
    elif version == 4:
        bc_new = bc_val[:]
        bc_new[n_D-1] /= 2
        k = np_rng.integers(1,6)
        load_config['loadstep'][0] = load_dict(n_D, bc_new, key)
        for m in np.arange(1,k+1):
            load_config['loadstep'].append(load_dict(n_D, ['x', 'x', 'x'], key))
        load_config['loadstep'].append(load_config['loadstep'][0] if key[-4:] == '_dot' else
                                       load_dict(n_D, bc_val, key))
    elif version == 5:
        k = np_rng.integers(1,5)
        bc_new = bc_val[:]
        bc_new[n_D-1] = -bc_new[n_D-1] if key[-4:] == '_dot' else 0.0
        load_config['loadstep'][0] = load_dict(n_D, bc_val, key, N = np_rng.integers(1,10))
        load_config['loadstep'].append(load_dict(n_D, bc_new, key, N = np_rng.integers(1,10)))
        for m in np.arange(1,k+1):
            load_config['loadstep'].append(load_dict(n_D, bc_val, key, N = np_rng.integers(1,10)))
            load_config['loadstep'].append(load_dict(n_D, bc_new, key, N = np_rng.integers(1,10)))

    return load_config


def material_setup(mat_config, np_rng):
    """
    Randomized material properties.

    The stiffness matrix needs to be positive definite.
    Here we ensure also a reasonable condition number because
    agreement with the analytical solution is better in
    that case.
    """
    unstable = True
    while unstable:
        C_11 = C_33 = 70.e9 + 20.e9 * np_rng.random()
        C_12 = C_13 = 60.e9 + 10.e9 * np_rng.random()
        C_44 = C_66 = 0.5 * (C_11 - C_12)

        C = np.array([[C_11, C_12, C_13,  0.0,  0.0,   0.0],
                      [C_12, C_11, C_12,  0.0,  0.0,   0.0],
                      [C_13, C_12, C_33,  0.0,  0.0,   0.0],
                      [ 0.0,  0.0,  0.0, C_44,  0.0,   0.0],
                      [ 0.0,  0.0,  0.0,  0.0,  C_44,  0.0],
                      [ 0.0,  0.0,  0.0,  0.0,  0.0,  C_66]])

        eigvals = np.linalg.eigvalsh(C)
        unstable = eigvals[0]/eigvals[-1] < 0.01
    print('eigenvalues', eigvals)

    # Set up material file
    mat_config['phase']['phase_0']['mechanical']['elastic'] = \
        {'type':'Hooke', 'C_11': C_11, 'C_12': C_12, 'C_44': C_44,
                         'C_33': C_33, 'C_13': C_13, 'C_66': C_66}
    mat_config = mat_config.material_add(phase = ['phase_0'],
                                         O = damask.Rotation.from_random(15,rng_seed=np_rng),
                                         homogenization = 'SX')
    return mat_config


@pytest.mark.parametrize('n_D', [2, 3])
@pytest.mark.parametrize('BC_type', ['u'])
def test_mesh_BC_coherence(res_path, copy_files, tmp_path, np_rng,
                           n_D, BC_type, petsc_version):
    copy_files(res_path, tmp_path, [f'mesh_{n_D}D.msh', f'load_{n_D}D.yaml', \
                                    'numerics.yaml'])

    load_config = damask.YAML.load(res_path/f'load_{n_D}D.yaml')
    mat_config = damask.ConfigMaterial.load(res_path/'material_base.yaml')

    bc_val = bc_setup(np_rng, n_D, BC_type)
    material_setup(mat_config, np_rng).save(tmp_path/'material.yaml')

    """ Reference solution (version 0, f_dot)"""
    BC_key = f'{BC_type}_dot'
    load_setup(load_config, np_rng, bc_val, n_D, BC_key, 0).save(tmp_path/f'load_v0_{n_D}D_{BC_key}.yaml')
    damask.util.run(f'damask_mesh -m material.yaml -l load_v0_{n_D}D_{BC_key}.yaml ' +
                    f'-g mesh_{n_D}D.msh -n numerics.yaml -j v0_{n_D}D_{BC_key}', wd = tmp_path)
    u_ref = damask.Result(tmp_path/f'v0_{n_D}D_{BC_key}.hdf5').view(increments=-1).get('u_n')
    u_ref[np.where(np.abs(u_ref) < 1.0e-9)] = 0.0

    """ Numerical solution"""
    BC_key = BC_type + '_dot'
    for t in ['rate']:
        for v in np.arange(1,5+1):
            load_config = damask.YAML.load(tmp_path/f'load_v0_{n_D}D_{BC_type}_dot.yaml')
            load_setup(load_config, np_rng, bc_val, n_D, BC_key, v).save(tmp_path/f'load_v{v}_{n_D}D_{BC_key}.yaml')
            damask.util.run(f'damask_mesh -m material.yaml -l load_v{v}_{n_D}D_{BC_key}.yaml ' +
                            f'-g mesh_{n_D}D.msh -n numerics.yaml -j v{v}_{n_D}D_{BC_key}', wd = tmp_path)
            u_n = damask.Result(tmp_path/f'v{v}_{n_D}D_{BC_key}.hdf5').view(increments=-1).get('u_n')
            u_n[np.where(np.abs(u_n) < 1.0e-12)] = 0.0
            assert(np.allclose(u_ref, u_n, rtol = 5.0e-6, atol = 1.0e-7) if v < 5 else
                   np.allclose(0.0, u_n, rtol = 1.0e-9, atol = 1.0e-12))
