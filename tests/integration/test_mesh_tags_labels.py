# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np
import damask

""" ----------------------------------------------------------------------------

    *** TEST FOR TAGS AND LABELS EQUIVALENCE ***

    Ensure the results of a simulation are the same whether boundary conditions
    are applied using the corresponding tags or labels.

    Parametrization:
     - Mesh file
     - Label to apply boundary condition
        - Use the label 'fixed' to fix degrees of freedom

    Randomization:
     - Orientations (in the material.yaml file)
     - Magnitude of the displacement boundary condition u_dot

    Resources (in resources/mesh_tags_labels):

    - load.yaml: load file (-l)
    - material_base.yaml: material file (-m; without orientations)
    - List of meshes (from all other mesh tests)
    - numerics.yaml: numerics configuration (-n)

    ------------------------------------------------------------------------ """

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'mesh_tags_labels'

def bc_setup(tmp_path, load_file, np_rng, n_D):
    """Randomized displacement boundary conditions."""
    # Randomly set displacement BC (v in 2D, w in 3D)
    dot_u = ['x', 'x', 'x']
    dot_u[int(n_D)-1] = np_rng.uniform(-1.e-2,+1.e-2)

    return dot_u


def load_setup(tmp_path, load_file, mesh, np_rng, dot_u, n_D, key, label):
    """Randomized displacement boundary conditions."""
    # Set up load file
    # Use key 'label' or 'tag'. If the latter, retrieve the tag corresponding
    # to 'label' from the mesh file
    load_config = damask.YAML.load(tmp_path/f'{load_file}.yaml')

    if key == 'by_label':
        load_config['loadstep'][0]['boundary_conditions']['mechanical'][1]['label'] = label
        load_config['loadstep'][0]['boundary_conditions']['mechanical'][1]['dot_u'] = dot_u
    else:
        tag = get_tag_from_label(tmp_path, mesh, label)
        rename_key_nested(load_config, 'label', 'tag')
        load_config['loadstep'][0]['boundary_conditions']['mechanical'][1]['tag']   = label
        load_config['loadstep'][0]['boundary_conditions']['mechanical'][1]['dot_u'] = dot_u

    load_config.save(tmp_path/f'{load_file}.yaml')


def material_setup(tmp_path, mat_file, np_rng):
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
    mat_config = damask.YAML.load(tmp_path/f'{mat_file}_base.yaml')
    mat_config['phase']['phase_0']['mechanical']['elastic'] = \
        {'type':'Hooke', 'C_11': C_11, 'C_12': C_12, 'C_44': C_44,
                         'C_33': C_33, 'C_13': C_13, 'C_66': C_66}
    mat_config.save(tmp_path/f'{mat_file}.yaml')
    mat_config = damask.ConfigMaterial.load(tmp_path/f'{mat_file}.yaml')
    mat_config = mat_config.material_add(phase = ['phase_0'],
                                         O = damask.Rotation.from_random(15,rng_seed=np_rng),
                                         homogenization = 'SX')
    mat_config.save(tmp_path/f'{mat_file}.yaml')


def damask_results(tmp_path, job):
    """Retrieve numerical solution from DADF5 file"""
    res = damask.Result(tmp_path/f'{job}.hdf5').view(increments=1)

    return res.get('u_n')


def get_tag_from_label(tmp_path, mesh, label):
    """ Find number of physical groups."""
    with open(tmp_path/f'{mesh}.msh', 'r') as f:
        while True:
            if f.readline().rstrip() == '$PhysicalNames':
                break

        n_PG = int(f.readline()) # number of tags/labels
        while n_PG > 0:
            l = f.readline().split()
            if l[2].replace('"','') == label:
                return int(l[1])
            n_PG -= 1

        raise ValueError(f'"{label}" not found')


def rename_key_nested(dictionary, old_key, new_key):
    """ Rename a key in a nested dictionary."""
    # https://www.geeksforgeeks.org/python/how-to-change-the-name-of-a-key-in-dictionary/
    for key in list(dictionary.keys()):
        if isinstance(dictionary[key], dict):
            rename_key_nested(dictionary[key], old_key, new_key)
        if key == old_key:
            dictionary[new_key] = dictionary.pop(old_key)


try:
    mark_old_PETSc = pytest.mark.xfail(pytest.petsc_version()<'3.24.0', reason='Not implemented')
except EnvironmentError:
    mark_old_PETSc = pytest.mark.xfail(False, reason='Assume working')
@pytest.mark.parametrize('n_D,   mesh_file,       label',
                         [('2', 'patch_tri',     'node_16'),
                          ('2', 'patch_tri',     'node_17'),
                          ('2', 'patch_quad',    'node_16'),
                          ('2', 'patch_quad',    'node_17'),
                          ('3', 'patch_tet',     'node_18'),
                          ('3', 'patch_tet',     'node_19'),
                          ('3', 'patch_hex',     'node_18'),
                          ('3', 'patch_hex',     'node_19'),
                          ('2', 'bicrystal_2D',  'tensionX'),
                          ('3', 'bicrystal_3D',  'tensionX'),
                          ('2', 'singleCrystal', 'tensionZ')])
def test_mesh_tags_labels(res_path, copy_files, tmp_path, np_rng,
                          n_D, mesh_file, label, petsc_version):
    copy_files(res_path, tmp_path)
    mesh = f'{mesh_file}'
    load = 'load'
    mat  = 'material'

    dot_u = bc_setup(tmp_path, load, np_rng, n_D)
    material_setup(tmp_path, mat, np_rng)

    for key in ['by_label', 'by_tag']:
        job = f'{mesh}_{key}'

        load_setup(tmp_path, load, mesh, np_rng, dot_u, n_D, key, label)

        """ Numerical solution """
        try:
            if n_D == 3 and dot_u.count('x') == 0 and petsc_version() < '3.24.1':
                pytest.xfail('solved in https://gitlab.com/petsc/petsc/-/merge_requests/8188')
        except EnvironmentError:
            pass

        damask.util.run(f'damask_mesh -m {mat}.yaml -l {load}.yaml -g {mesh}.msh ' +
                        f'-n numerics.yaml -j {job}', wd = tmp_path)

    # Results comparison
    with_label = damask_results(tmp_path, f'{mesh}_by_label'),
    with_tag   = damask_results(tmp_path, f'{mesh}_by_tag'),
    assert(np.allclose(with_label, with_tag, rtol = 1.0e-14, atol = 1.0e-14))
