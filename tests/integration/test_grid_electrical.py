# SPDX-License-Identifier: AGPL-3.0-or-later

"""Prototype DAMASK electrical sphere regression-test harness.

This script is meant as a clean discussion/prototype artifact for a DAMASK test:

1. create a cubic grid with a centered spherical inclusion,
2. write minimal material and load YAML files,
3. write the analytical reference electric field,
4. optionally run DAMASK_grid,
5. optionally compare a DAMASK electric-field export against the reference.

"""
import subprocess
from pathlib import Path

import numpy as np
import damask
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


LOAD_YAML = """\
solver: {{mechanical: spectral_polarization, thermal: spectral, electrical: spectral}}
loadstep:
- discretization: {{t: {time}, N: 1}}
  boundary_conditions:
    mechanical:
      F:
      - [x, 0.0, 0.0]
      - [0.0, x, 0.0]
      - [0.0, 0.0, x]
      P:
      - [0.0, x, x]
      - [x, 0.0, x]
      - [x, x, 0.0]
    electrical:
      E: [{e0}, 0.0, 0.0]
      J: [x, x, x]
"""


def material_yaml(sigma_matrix: float, sigma_inclusion: float) -> str:
    return f"""\
homogenization:
  single:
    mechanical:
      type: pass
    thermal:
      type: pass
      output: [T]
    electrical:
      type: pass
      output: [E, J]
    N_constituents: 1
phase:
  matrix:
    lattice: cF
    rho: 1740.0
    mechanical:
      output: [O]
      elastic: {{type: Hooke, C_11: 50e9, C_12: 20e9, C_44: 15e9}}
      eigen:
      - {{type: thermalexpansion, Alpha_11: 25.0e-6, Alpha_33: 26.0e-6}}
    thermal:
      C_p: 1051.7
      K_11: 160.0
      K_33: 156.0
      source: [{{type: Joule}}]
    electrical:
      sigma: {sigma_matrix:.8e}
      type: iso
      output: [E, J]
  inclusion:
    lattice: cF
    rho: 1740.0
    mechanical:
      output: [O]
      elastic: {{type: Hooke, C_11: 50e9, C_12: 20e9, C_44: 15e9}}
      eigen:
      - {{type: thermalexpansion, Alpha_11: 25.0e-6, Alpha_33: 26.0e-6}}
    thermal:
      C_p: 1051.7
      K_11: 160.0
      K_33: 156.0
      source: [{{type: Joule}}]
    electrical:
      sigma: {sigma_inclusion:.8e}
      type: iso
      output: [E, J]
material:
- constituents:
  - phase: matrix
    O: [1.0, 0.0, 0.0, 0.0]
    v: 1.0
  homogenization: single
- constituents:
  - phase: inclusion
    O: [1.0, 0.0, 0.0, 0.0]
    v: 1.0
  homogenization: single
"""


def analytical_electric_field(cells,size,radius,sigma_inclusion,sigma_matrix,e0):
    coordinates = damask.grid_filters.coordinates0_point(cells,size,origin=-size/2.)
    r = np.linalg.norm(coordinates,axis=-1)
    r_safe = np.where(r == 0.0, 1.0, r)
    dipole_factor = (sigma_inclusion - sigma_matrix) / (sigma_inclusion + 2.0 * sigma_matrix)
    inclusion_factor = 3.0 * sigma_matrix / (sigma_inclusion + 2.0 * sigma_matrix)

    field = np.zeros(tuple(cells)+(3,), dtype=float)
    inside = r <= radius
    outside = ~inside
    field[inside, 0] = inclusion_factor * e0

    if np.any(outside):
        ro = r_safe[outside]
        xo = coordinates[...,0][outside]
        yo = coordinates[...,1][outside]
        zo = coordinates[...,2][outside]
        radius_over_r_cubed = (radius / ro) ** 3
        r2 = ro**2
        field[outside, 0] = e0 * (
            1.0
            + dipole_factor
            * radius_over_r_cubed
            * (3.0 * xo**2 / r2 - 1.0)
        )
        field[outside, 1] = (e0 * dipole_factor * radius_over_r_cubed * 3.0 * xo * yo / r2)
        field[outside, 2] = (e0 * dipole_factor * radius_over_r_cubed * 3.0 * xo * zo / r2)

    center = tuple(cells//2)
    if r[center] == 0.0:
        field[center] = [inclusion_factor * e0, 0.0, 0.0]
    return field


def create_geometry(path, cells, size, radius):
    geom = damask.GeomGrid(material=np.zeros(cells, dtype=int), size=size). \
           add_primitive(dimension=2.0 * radius, center = 0.5*size,exponent=1)
    geom.initial_conditions['T'] = 300.0
    geom.save(path)


def compare_fields(numerical, reference):
    diff = numerical - reference
    rel_l2 = np.linalg.norm(diff.ravel()) / np.linalg.norm(reference.ravel())
    max_abs = float(np.max(np.abs(diff)))
    return float(rel_l2), max_abs


def centerline(field: np.ndarray) -> np.ndarray:
    """Return an x-line through the grid center.

    For even grids the exact symmetry line lies between cells. Average the four
    nearest y/z cell-center lines so transverse components respect symmetry
    better than a single off-center line.
    """
    n = field.shape[0]
    lo = n // 2 - 1
    hi = n // 2
    return 0.25 * (
        field[:, lo, lo, :]
        + field[:, lo, hi, :]
        + field[:, hi, lo, :]
        + field[:, hi, hi, :]
    )


def make_plot(numerical, reference, length, radius, out: Path):
    n = reference.shape[0]
    dx = length / n
    coords = (np.arange(n, dtype=float) + 0.5) * dx - 0.5 * length
    reference_line = centerline(reference)
    numerical_line = centerline(numerical)

    components = [(0, "E1"), (1, "E2 / E3")]
    fig, axes = plt.subplots(len(components), 1, figsize=(8, 6), sharex=True)
    for axis, (component, label) in zip(axes, components):
        axis.plot(coords, reference_line[:, component], "k-", label="analytical")
        axis.plot(coords, numerical_line[:, component], "r--", label="DAMASK")
        axis.axvline(-radius, color="0.5", ls=":")
        axis.axvline(radius, color="0.5", ls=":")
        axis.set_ylabel(label)
        axis.grid(True, alpha=0.3)
    axes[-1].set_xlabel("x - center")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def test_analytical(tmp_path):
    geom = 'inclusion.vti'
    material = 'material.yaml'
    load = 'load.yaml'
    job = 'job'

    cells = np.array([16]*3)
    length=1.0
    size = np.array([length]*3)
    radius_over_length=1.0 / 16.0
    sigma_matrix=3.0e7
    sigma_inclusion=2.0e7
    e0=1.0
    time=0.002
    radius = length * radius_over_length

    create_geometry(tmp_path/geom, cells, size, radius)
    m = damask.ConfigMaterial(material_yaml(sigma_matrix, sigma_inclusion))
    m.save(tmp_path/material)
    g = damask.LoadcaseGrid(LOAD_YAML.format(time=time, e0=e0))
    g.save(tmp_path/load)
    damask.util.run(f'damask_grid -l {load} -g {geom} -m {material} -j {job}',wd=tmp_path)

    analytical = analytical_electric_field(cells,size,radius,sigma_inclusion,sigma_matrix,e0)
    numerical = damask.Result(tmp_path/f'{job}.hdf5').view(increments=-1,phases=False).place('E')
    numerical = damask.grid_filters.unravel(numerical,cells)

    rel_l2, max_abs = compare_fields(numerical, analytical)
    make_plot(numerical, analytical, length, radius,
              tmp_path/'electrostatic_sphere_line_compare.png')
    assert rel_l2 < 1.e-2
    assert max_abs < 1.e-1
