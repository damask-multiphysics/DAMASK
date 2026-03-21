# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np
import h5py
from vtkmodules.vtkFiltersVerdict import vtkCellSizeFilter

import damask


""" ----------------------------------------------------------------------------

    *** VTK HIGH ORDER ELEMENT CONFORMANCE ***

    Check that the geometry data in the output HDF5 (node coordinates 'x_n' and
    connectivity matrix 'T_e') provide the correct values in order to build a
    conforming Lagrange VTK element.

    Parametrization:
     - Element type (polytope):
       - Triangle (TRIA)
       - Quadrilateral (QUAD)
       - Tetrahedron (TETR)
       - Hexahedron (HEXA)
     - FE approximation order
       - 1 to 5

    Meshes:
     - Single element in [-1,+1]^d
     - Multiple elements in [-1,+1]^d (single element subdivided)

    Exact values for the number of nodes and the area/volume of the elements are
    also provided.

    Resources (in resources/mesh_high_order_vtk):
    - Load (-l)
      - load.yaml
    - Material (-m)
      - material.yaml
    - Mesh (-g)
      - {TRIA/QUAD/TETR/HEXA}1.msh (single element)
      - {TRIA/QUAD/TETR/HEXA}m.msh (multiple element)
    - Numerics configuration (-n)
      - numerics.yaml

    ------------------------------------------------------------------------ """


"""Known element data."""
nodes = {'TRIA': [ 3,  6, 10,  15,  21],                                                            # Number of nodes per order
         'QUAD': [ 4,  9, 16,  25,  36],
         'TETR': [ 4, 10, 20,  35,  56],
         'HEXA': [ 8, 27, 64, 125, 216]}

sizes = {'TRIA': 2,                                                                                 # Area / volume
         'QUAD': 4,
         'TETR': 4/3,
         'HEXA': 8}


@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'mesh_high_order_vtk'


@pytest.mark.parametrize('polytope',  list(nodes.keys()))
@pytest.mark.parametrize('order', [1, 2, 3, 4, 5])                                                  # FE approximation space order
@pytest.mark.parametrize('elements', ['single', 'multi'])                                           # Single or multi-element mesh
def test_mesh_high_order_VTK(res_path, copy_files, tmp_path, polytope, order, elements):
    mesh = f'{polytope}{"1" if elements == "single" else "m"}.msh'
    job  = f'{mesh}_o{order}'
    copy_files(res_path, tmp_path, ['material.yaml','load.yaml',mesh])

    """Numerical solution."""
    numerics = damask.YAML.load(res_path/'numerics.yaml')
    numerics['solver']['mesh']['p_s'] = order
    numerics['solver']['mesh']['p_i'] = order
    numerics.save(tmp_path/'numerics.yaml')

    damask.util.run(f'damask_mesh -l load.yaml -g {mesh} -m material.yaml ' +
                    f'-n numerics.yaml -j {job}', wd = tmp_path)

    """Create unstructured grid data."""
    with h5py.File(tmp_path/f'{job}.hdf5', 'r') as f:

        PETSc_version = f'{f.attrs["PETSc_version_major"]}.{f.attrs["PETSc_version_minor"]}.{f.attrs["PETSc_version_subminor"]}'
        if damask.util.version(PETSc_version) < '3.24.1': pytest.xfail()

        x_n = f['/geometry']['x_n'][()]                                                             # nodal coordinates
        T_e = f['/geometry']['T_e'][()]                                                             # element connectivity

    assert T_e.shape[1] == nodes[polytope][int(order)-1]                                            # nodes per cell

    vtk = damask.VTK.from_unstructured_grid(x_n,T_e,polytope)
    unstructured_grid = vtk.vtk_data

    sizeFilter = vtkCellSizeFilter()
    sizeFilter.SetInputData(unstructured_grid)
    sizeFilter.Update()

    geom = 'Area' if unstructured_grid.GetCell(0).GetCellDimension() == 2 else 'Volume'
    vals = sizeFilter.GetOutput().GetCellData().GetArray(geom)
    size = sum([vals.GetValue(i) for i in range(vals.GetNumberOfTuples())])
    assert np.isclose(size, sizes[polytope])
