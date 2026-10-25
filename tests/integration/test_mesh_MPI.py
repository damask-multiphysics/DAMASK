# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np

import damask

""" ----------------------------------------------------------------------------

    *** COMPARE RESULTS OF PARALLEL SIMULATIONS ***

    Ensure the reults of a simulation are the same for any number of MPI ranks.

    Parametrization:
     - Mesh file element and fineness
     - Element order (1/2/3)

    Randomization:
     - Orientations (in the material.yaml file)
     - Magnitude of the displacement boundary condition

    Resources (in resources/mesh_MPI):

    - load_[tq/colZ].yaml: load file (-l)
    - material_base.yaml: material file (-m; without orientations)
    - meshes:
        - tri/quad   (triangles/quadrilaterals)
        - colZ[5/8g] (tetrahedra/hexahedra)
    - numerics.yaml: numerics configuration (-n)

    ------------------------------------------------------------------------ """


@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'mesh_MPI'


def numerics_setup(res_path,tmp_path,p_s,simplex):
    """Update shape function and integration orders in numerics file"""
    numerics_config = damask.YAML.load(res_path/'numerics.yaml')
    numerics_config['solver']['mesh']['p_s'] = p_s
    numerics_config['solver']['mesh']['p_i'] = np.maximum(2*(p_s-1),1) if simplex else \
                                               p_s+1
    numerics_config.save(tmp_path/'numerics.yaml')


def load_setup(load_config,np_rng,colZ):
    """Randomized boundary conditions."""
    val = np_rng.choice([-1,+1])*np_rng.uniform(1,5)*10**(-3)

    load_config['loadstep'][0]['boundary_condition']['mechanical'].append(
        {'label': 'surf_top'      if colZ else 'e_2',
         'u':     ['x', 'x', val] if colZ else [val, 'x', 'x']})

    return load_config


def material_setup(mat_config,np_rng):
    """
    Randomized material properties.

    The stiffness matrix needs to be positive definite.
    Here we ensure also a reasonable condition number because
    agreement with the analytical solution is better in
    that case.
    """
    unstable = True
    while unstable:
        C_11 = 80.e9 + 25.e9 * np_rng.random()
        C_12 = 60.e9 + 15.e9 * np_rng.random()
        C_44 = 35.e9 + 10.e9 * np_rng.random()
        C_33 = 70.e9 + 15.e9 * np_rng.random()
        C_13 = 60.e9 + 10.e9 * np_rng.random()
        C_66 = 0.5 * (C_11 - C_12)

        C = np.array([[C_11, C_12, C_13,  0.0,  0.0,   0.0],
                      [C_12, C_11, C_13,  0.0,  0.0,   0.0],
                      [C_13, C_13, C_33,  0.0,  0.0,   0.0],
                      [ 0.0,  0.0,  0.0, C_44,  0.0,   0.0],
                      [ 0.0,  0.0,  0.0,  0.0,  C_44,  0.0],
                      [ 0.0,  0.0,  0.0,  0.0,  0.0,  C_66]])

        eigvals = np.linalg.eigvalsh(C)
        unstable = eigvals[0]/eigvals[-1] < 0.01

    # Set up material file
    mat_config['phase']['phase_0']['mechanical']['elastic'] = \
        {'type':'Hooke', 'C_11': C_11, 'C_12': C_12, 'C_44': C_44,
                         'C_33': C_33, 'C_13': C_13, 'C_66': C_66}
    mat_config = mat_config.material_add(phase = ['phase_0'],
                                         O = damask.Rotation.from_random(8,rng=np_rng),
                                         homogenization = 'SX')
    return mat_config


def ccc(x, y):
    """Concordance correlation coefficient."""
    x_bar, y_bar = np.mean(x), np.mean(y)
    s_xx, s_yy = np.var(x), np.var(y)
    s_xy = np.mean((x - x_bar) * (y - y_bar))
    return 2.0 * s_xy / (s_xx + s_yy + (x_bar - y_bar)**2)


try:
    mark_old_PETSc = pytest.mark.xfail(pytest.petsc_version()<'3.23.0', reason='Not working')
except EnvironmentError:
    mark_old_PETSc = pytest.mark.xfail(False, reason='Assume working')


@pytest.mark.parametrize('mesh,    p_s',
                         [('tri',    1),   # TRI
                          ('tri',    2),   # TRI
                          ('tri',    3),   # TRI
                          ('quad',   1),   # QUAD
                          ('quad',   2),   # QUAD
                          ('quad',   3),   # QUAD
                          ('colZ5g', 1),   # TET
                          ('colZ5g', 2),   # TET
                          ('colZ8g', 1),   # HEX
                          ('colZ8g', 2)])  # HEX
def test_mesh_MPI(res_path,copy_files,tmp_path,mpi_launcher,assert_allclose,np_rng,mesh,p_s):
    """ Check size/dimension and values of output data."""

    copy_files(res_path,tmp_path,[f'{mesh}.msh'])

    load = 'load_colZ' if mesh.endswith('g') else 'load_tq'
    material = 'material'

    simplex = mesh.startswith(('tri','colZ16g'))
    load_setup(damask.YAML.load(res_path/f'{load}.yaml'),np_rng,
                                load.endswith('Z')).save(tmp_path/f'{load}.yaml')
    material_setup(damask.ConfigMaterial.load(res_path/'material_base.yaml'),
                   np_rng).save(tmp_path/'material.yaml')

    numerics_setup(res_path,tmp_path,p_s,simplex)

    # Reference (serial)
    job = f'n_1_{mesh}_o{p_s}'
    damask.util.run(f'{mpi_launcher} -n 1 ' +
                    f'damask_mesh -l {load}.yaml -g {mesh}.msh -m {material}.yaml ' +
                    f'-n numerics.yaml -j {job}', wd = tmp_path)
    ref_r = damask.Result(tmp_path/f'{job}.hdf5').view(increments=-1)

    # Geometry data
    ref_x0_n, ref_x0_p, ref_T_e = ref_r.coordinates0_node, ref_r.coordinates0_point, ref_r.connectivity
    ref_shape  = np.array([np.shape(ref_x0_n), np.shape(ref_x0_p), np.shape(ref_T_e)])

    # Sorting indices (nodes & IPs)
    ref_idx_n = np.lexsort(tuple(ref_x0_n[:, i] for i in range(ref_x0_n.shape[1])))
    ref_idx_p = np.lexsort(tuple(ref_x0_p[:, i] for i in range(ref_x0_p.shape[1])))

    # Displacements, deformation gradient & stresses
    ref_u, ref_F, ref_P = ref_r.get('u_n'), ref_r.get('F'), ref_r.get('P')

    # Parallel
    for rank in [2,3,4]:
        job = f'n_{rank}_{mesh}_o{p_s}'
        damask.util.run(f'{mpi_launcher} -n {rank} ' +
                        f'damask_mesh -l {load}.yaml -g {mesh}.msh -m {material}.yaml ' +
                        f'-n numerics.yaml -j {job}', wd = tmp_path)
        r = damask.Result(tmp_path/f'{job}.hdf5').view(increments=-1)

        # Geometry data
        x0_n, x0_p, T_e = r.coordinates0_node, r.coordinates0_point, r.connectivity
        n_shape = np.array([np.shape(x0_n), np.shape(x0_p), np.shape(T_e)])

        # Assert correct ouptput shape (no ghost/repeated values)
        assert np.all(ref_shape == n_shape)

        # Sorting indices
        idx_n = np.lexsort(tuple(x0_n[:, i] for i in range(x0_n.shape[1])))
        idx_p = np.lexsort(tuple(x0_p[:, i] for i in range(x0_p.shape[1])))

        # Displacements, deformation gradient & stresses
        u, F, P = r.get('u_n'), r.get('F'), r.get('P')

        # Assert geometry data
        assert_allclose(ref_x0_n[ref_idx_n],x0_n[idx_n])
        assert_allclose(ref_x0_p[ref_idx_p],x0_p[idx_p])

        # Assert solution values
        rho_c_u = ccc(ref_u[ref_idx_n].flatten(),u[idx_n].flatten())
        rho_c_F = ccc(ref_F[ref_idx_p].flatten(),F[idx_p].flatten())
        rho_c_P = ccc(ref_P[ref_idx_p].flatten(),P[idx_p].flatten())
        assert (np.array([rho_c_F, rho_c_u, rho_c_P]) > 0.99).all()
