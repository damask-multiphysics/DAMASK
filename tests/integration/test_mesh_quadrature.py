# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np

import damask

""" ----------------------------------------------------------------------------

    *** QUADRATURE ORDER EQUIVALENCE ***

    Ensure simulation results are correct (against an analytical solution) when
    using different combinations of element and integration order, including
    reduced integration if appropriate.

    - test_mesh_quadrature
        Same mesh (#cells), increasing element and quadrature orders.
    - test_mesh_quadrature_ho (high order)
        With higher order elements, the same results can be obtained with
        coarser meshes.
        The parametrized variable 'size' indicates the number of nodes along
        the arch edges.

    Parametrization:
     - Element type:
       - Triangle (tri)
       - Quadrilateral (quad)
     - FE approximation order
       - 1 to 5 (tri & quad)
     - Integration/quadrature order
       - Full
       - Reduced
     - Boundary condition
       - Magnitude of the displacement
     - Elastic constants
     - Mesh coarseness

    Geometry:
     - One quarter, curved cantilever beam with inner/outer radii of 5/10 m.

    - Meshes:
       - test_mesh_quadrature:
            760 elements (tri & quad) / 20 x 38 (quads)
       - test_mesh_quadrature_ho (quads)
            order 2: 6 x 14
            order 3: 6 x 12
            order 4: 6 x 10
            order 5: 4 x 8

    Reference solution:
     - Theory of Elasticity.
       S. Timoshenko & J. N. Goodier
       McGraw-Hill, second Edition, 1951.
       Problem #31. Bending of a Curved Bar by a Force at the End.

     - The Finite Element Method: Its Basis and Fundamentals.
       O. C. Zienkiewicz, R. L. Taylor & J. Z. Zhu
       Butterworth-Heinemann, seventh edition, 2013.
       Example 7.6. Circular beam subjected to end shear

    Resources (in resources/mesh_quadrature):
    - Load (-l)
      - load_curved.yaml
    - Material (-m)
      - material_base.yaml (without elastic constants)
    - Mesh (-g)
      - curved_{tri/quad}.msh
    - Numerics configuration (-n)
      - numerics.yaml

    Runtime:
    - test_mesh_quadrature:     ~ 1 min 15 s
    - test_mesh_quadrature_ho:  < 10 s
    ------------------------------------------------------------------------ """


@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'mesh_quadrature'


def load_setup(load_config, tmp_path, np_rng):
    """Randomized displacement boundary conditions."""
    u_BC = -1.0e-3*np_rng.uniform(0.5,5.0)
    load_config['loadstep'][0]['boundary_condition']['mechanical'].append(
            {'label': 'bot', 'u': [u_BC, 'x', 'x']})
    load_config.save(tmp_path/'load_curved.yaml')

    return u_BC


def material_setup(file_in, file_out, np_rng):
    """
    Randomized material properties.

    Use an isotropic material.

    E  = (C_11-C_12)*(C_11+2*C_12)/(C_11+C_12)
    nu = C_12/(C_11+C_12)

    C_11 = E*(1-nu)/((1+nu)*(1-2*nu))
    C_12 = E*nu/((1+nu)*(1-2*nu))
    C_44 = nu/(2*(1+nu))

    The stiffness matrix needs to be positive definite.
    We additionally constrain its condition number
    to max(eigvals)/min(eigvals) < 100,
    since well-conditioned systems typically exhibit
    closer agreement with the analytical solution.
    """
    unstable = True
    while unstable:
        E    = np_rng.uniform(0.5,10)
        nu   = np_rng.uniform(0,0.375)
        C_11 = E*(1-nu)/((1+nu)*(1-2*nu))
        C_12 = E*nu/((1+nu)*(1-2*nu))
        C_44 = E/(2*(1+nu))

        C = np.array([[C_11, C_12, C_12,  0.0,  0.0,   0.0],
                      [C_12, C_11, C_12,  0.0,  0.0,   0.0],
                      [C_12, C_12, C_11,  0.0,  0.0,   0.0],
                      [ 0.0,  0.0,  0.0, C_44,  0.0,   0.0],
                      [ 0.0,  0.0,  0.0,  0.0,  C_44,  0.0],
                      [ 0.0,  0.0,  0.0,  0.0,  0.0,  C_44]])
        eigvals = np.linalg.eigvalsh(C)
        unstable = np.min(eigvals)/np.max(eigvals) < 0.01

    # Set up material file
    mat_config = damask.ConfigMaterial.load(file_in)
    mat_config['phase']['phase_0']['mechanical']['elastic'] = \
        {'type':'Hooke', 'C_11': C_11, 'C_12': C_12, 'C_44': C_44}
    mat_config.save(file_out)

    return np.array([E, nu])


def numerics_setup(res_path, tmp_path, p_s, r_int):
    """Update shape function and integration orders in numerics file."""
    numerics_config = damask.YAML.load(res_path/'numerics.yaml')
    numerics_config['solver']['mesh']['p_s'] = p_s
    if r_int: numerics_config['solver']['mesh']['reduced_integration'] = r_int
    numerics_config.save(tmp_path/'numerics.yaml')


""" Node indices in the mesh """
def nodes_idx(el, order):
    """Indices of nodes belonging to each boundary.

    The array 'st' containst the nodes corresponding to the four corner vertices

        (a,0), (b,0), (0,a), (0,b),         a = 5, b = 10

    The nodes along the four boundaries 'inn', 'out', 'bot', 'top' are then numbered
    in a contiguous way after these four.

    The mesh is the same for all simulations, and so is the number of nodes in each
    boundary (array 'nn'). Higher-order nodes are not considered (it would simply add
    more nodes for the assertion of results).

        st[0]: (a, 0)
        ...
        st[0] + 4: first node within 'inn' boundary
        ...
        st[0] + 4 + nn[0]: last+1 node within 'inn' boundary / first node within 'out' boundary

    "within" not considering the corner nodes.

    """
    if el == 'tri':
        st = np.array([0, 0, 760, 2280, 4560])
        nn = np.array([17, 35, 11, 11])
    elif el == 'quad':
        st = np.array([0, 760, 3040, 6840, 12160])
        nn = np.array([37, 37, 19, 19])
    else:
        st = np.array([84, 288, 540, 512])
        nn = np.array([[13, 13, 5, 5],      # order 2
                       [11, 11, 5, 5],
                       [ 9,  9, 5, 5],
                       [ 7,  7, 3, 3]])     # order 5

    ss = order - 2 if el == 'ho' else order - 1
    cs = st[ss] + 4 + np.concatenate(([0], np.cumsum(nn[ss] if el == 'ho' else nn)))
    idx = {
        'inn': np.concatenate(([st[ss]+0], np.arange(cs[0], cs[1]), [st[ss]+2])),
        'out': np.concatenate(([st[ss]+1], np.arange(cs[1], cs[2]), [st[ss]+3])),
        'bot': np.concatenate(([st[ss]+0], np.arange(cs[2], cs[3]), [st[ss]+1])),
        'top': np.concatenate(([st[ss]+2], np.arange(cs[3], cs[4]), [st[ss]+3])),
    }

    return idx


def displ_analytical(elast_prop, u0, edge, x_n, chop = 1.0e-12):
    """ Analytical solution."""
    aa, bb = 5.0, 10.0
    EE, nu = elast_prop[0], elast_prop[0]
    NN = aa**2 - bb**2 + (aa**2 + bb**2) * np.log(bb / aa)
    PP = -u0 * NN * EE / np.pi / (aa**2 + bb**2)
    PN = PP / (EE * NN)
    KK = PN * (1/2*(1-3*nu)*aa**2 - bb**2*(1+nu)/2 - (aa**2+bb**2)*(1-nu)*np.log(aa))

    rr, tt = r_t_analytical(x_n, edge, aa, bb)

    u_r = PN * (
            (1/2*(1-3*nu)*rr**2 - aa**2*bb**2*(1+nu)/(2*rr**2)
             - (aa**2+bb**2)*(1-nu)*np.log(rr)) * np.sin(tt)
            + (aa**2+bb**2) * (2*tt-np.pi) * np.cos(tt)
        ) - KK * np.sin(tt)

    u_t = - PN * (
            (1/2*(5+nu)*rr**2 - aa**2*bb**2*(1+nu)/2/rr**2
             + (aa**2+bb**2)*((1-nu)*np.log(rr)+(1+nu))) * np.cos(tt)
            + (aa**2+bb**2) * (2*tt-np.pi) * np.sin(tt)
        ) - KK * np.cos(tt)


    u_r[np.abs(u_r) < chop] = 0.0
    u_t[np.abs(u_t) < chop] = 0.0

    return u_r, u_t


def r_t_analytical(x_n, edge, a, b):
    """ Polar coordinates (r,theta)."""
    if edge == 'bot':
        r, t = x_n[:,0], 0.0
    elif edge == 'top':
        r, t = x_n[:,1], np.pi/2
    elif edge == 'inn':
        ang = x_n[:,1]/a
        ang[np.abs(ang) > 1.0] = 1.0
        r, t = a, np.asin(ang)
    elif edge == 'out':
        ang = x_n[:,1]/b
        ang[np.abs(ang) > 1.0] = 1.0
        r, t = b, np.asin(ang)

    return r, t


def displ_damask(edge, x_n, u_n, chop = 1.0e-12):
    """ DAMASK numerical solution."""
    a, b = 5.0, 10.0

    if edge == 'bot':
        u_r = u_n[:,0]
        u_t = u_n[:,1]
    elif edge == 'top':
        u_r = u_n[:,1]
        u_t = u_n[:,0]
    elif edge == 'inn':
        ang = x_n[:,1]/a
        ang[np.abs(ang) > 1.0] = 1.0
        ang = np.asin(ang)
        u_r = +np.cos(ang) * u_n[:,0] + np.sin(ang) * u_n[:,1]
        u_t = -np.sin(ang) * u_n[:,0] + np.cos(ang) * u_n[:,1]
    elif edge == 'out':
        ang = x_n[:,1]/b
        ang[np.abs(ang) > 1.0] = 1.0
        ang = np.asin(ang)
        u_r = +np.cos(ang) * u_n[:,0] + np.sin(ang) * u_n[:,1]
        u_t = -np.sin(ang) * u_n[:,0] + np.cos(ang) * u_n[:,1]

    u_r[np.abs(u_r) < chop] = 0.0
    u_t[np.abs(u_t) < chop] = 0.0

    return u_r, u_t


def damask_results(tmp_path, job):
    """Retrieve numerical solution from DADF5 file."""
    r = damask.Result(tmp_path/f'{job}.hdf5')
    x_n = r.coordinates0_node
    u_n = r.view(increments=-1).get('u_n')

    return x_n, u_n


@pytest.mark.parametrize(  'el,    p_s,  r_int',
                         [('tri',  i+1,  False)  for i in range(5)] +
                         [('quad',  1,   True)]                     +
                         [('quad', i+1,  False)  for i in range(5)])
>>>>>>> 4e9f1efc9 (Add test (II): high order elements)
def test_mesh_quadrature(res_path, copy_files, tmp_path, np_rng,
                         el, p_s, r_int):
    mesh = f'curved_{el}'
    load = 'load_curved'
    job  = f'{mesh}_o{p_s}' + ('r' if r_int else 'f')

    copy_files(res_path, tmp_path, [f'{mesh}.msh', 'numerics.yaml'])

    u_BC = load_setup(damask.YAML.load(res_path/f'{load}.yaml'), tmp_path, np_rng)
    elast_prop = material_setup(res_path/'material_base.yaml', tmp_path/'material.yaml', np_rng)
    numerics_setup(res_path, tmp_path, p_s, r_int)

    damask.util.run(f'damask_mesh -m material.yaml -l {load}.yaml -g {mesh}.msh ' +
                    f'-n numerics.yaml -j {job}', wd = tmp_path)
    x_n, u_n = damask_results(tmp_path, job)
    idx = nodes_idx(el, p_s)

    for edge in ['inn', 'out', 'bot', 'top']:
        u_r_D, u_t_D = displ_damask(edge, x_n[idx[edge]], u_n[idx[edge]])
        u_r,   u_t   = displ_analytical(elast_prop, u_BC, edge, x_n[idx[edge]])

        if el == 'quad' and p_s == 1 and r_int:                                                     # hourglassing
            # In this particular case, every other node has the correct displacement,
            # while the others are off
            u_r,   u_t   = u_r[::2],   u_t[::2]
            u_r_D, u_t_D = u_r_D[::2], u_t_D[::2]

        assert np.allclose(u_r_D, u_r, rtol = 1.0e-5, atol = 5.0e-3)
        assert np.allclose(u_t_D, u_t, rtol = 1.0e-5, atol = 5.0e-3)


@pytest.mark.parametrize(   'el,   p_s,   size',
                         [('quad', i+2,  15-2*i) for i in range(4)])
def test_mesh_quadrature_ho(res_path, copy_files, tmp_path, np_rng,
                            el, p_s, size):
    mesh = f'curved_{el}_ho_n{size}'
    load = 'load_curved'
    job  = f'curved_{el}_ho_{p_s}f'

    copy_files(res_path, tmp_path, [f'{mesh}.msh', 'numerics.yaml'])

    u_BC = load_setup(damask.YAML.load(res_path/f'{load}.yaml'), tmp_path, np_rng)
    elast_prop = material_setup(res_path/'material_base.yaml', tmp_path/'material.yaml', np_rng)
    numerics_setup(res_path, tmp_path, p_s, False)

    damask.util.run(f'damask_mesh -m material.yaml -l {load}.yaml -g {mesh}.msh ' +
                    f'-n numerics.yaml -j {job}', wd = tmp_path)
    x_n, u_n = damask_results(tmp_path, job)
    idx = nodes_idx('ho', p_s)

    for edge in ['inn', 'out', 'bot', 'top']:
        u_r_D, u_t_D = displ_damask(edge, x_n[idx[edge]], u_n[idx[edge]])
        u_r,   u_t   = displ_exact(elast_prop, u_BC, edge, x_n[idx[edge]])

        assert np.allclose(u_r_D, u_r, rtol = 1.0e-5, atol = 5.0e-3)
        assert np.allclose(u_t_D, u_t, rtol = 1.0e-5, atol = 5.0e-3)
