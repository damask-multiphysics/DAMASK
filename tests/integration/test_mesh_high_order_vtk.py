# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np
import h5py
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkFiltersVerdict import vtkCellSizeFilter
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid, vtkCellArray, \
                                          vtkLagrangeTriangle, vtkLagrangeTetra
from vtkmodules.util.numpy_support import numpy_to_vtk as np_to_vtk, \
                                          numpy_to_vtkIdTypeArray as np_to_vtkId

import damask


""" ----------------------------------------------------------------------------

    *** VTK HIGH ORDER ELEMENT CONFORMANCE ***

    Check that the geometry data in the output HDF5 (node coordinates 'x_n' and
    connectivity matrix 'T_e') provide the correct values in order to build a
    conforming Lagrange VTK element.

    Parametrization:
     - Element type (polytope):
       - Triangle (TRI)
       - Tetrahedron (TET)
     - FE approximation order
       - 1 to 5

    Meshes:
     - Single element in [-1,+1]^d
     - Multiple elements in [-1,+1]^d (single element subdivided)

    Exact values for the number of nodes and the area/volume of the elements are
    also provided.

    Resources (in resources/mesh_high_order_vtk):

        load.yaml       : load file (-l; no deformation)
        material.yaml   : material file (-m)
        TRI/TET.msh     : mesh file (-g)
        numerics.yaml   : numerics configuration (-n)

    ------------------------------------------------------------------------ """


"""Known element data."""
nodes = {'TRI': [ 3,  6, 10, 15, 21],                                                               # Number of nodes per order
         'TET': [ 4, 10, 20, 35, 56]}

sizes = {'TRI': 2,                                                                                  # Area / volume
         'TET': 4/3}


@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'mesh_high_order_vtk'


@pytest.mark.parametrize('polytope',  ['TRI', 'TET'])                                               # Polytope type
@pytest.mark.parametrize('order', [1, 2, 3, 4, 5])                                                  # FE approximation space order
@pytest.mark.parametrize('nElements', ['single', 'multi'])                                          # Single or multi-element mesh
def test_mesh_high_order_VTK(res_path, copy_files, tmp_path, polytope, order, nElements):
    copy_files(res_path, tmp_path)
    load = 'load'
    mesh = f'{polytope}1' if nElements == 'single' else f'{polytope}m'
    mat  = 'material'
    num  = 'numerics'
    job  = f'{mesh}_o{order}'

    """Numerical solution."""
    numerics = damask.YAML.load(tmp_path/f'{num}.yaml')
    numerics['solver']['mesh']['p_s'] = order
    numerics['solver']['mesh']['p_i'] = order
    numerics.save(tmp_path/f'{num}.yaml')

    damask.util.run(f'damask_mesh -l {load}.yaml -g {mesh}.msh -m {mat}.yaml ' +
                    f'-n {num}.yaml -j {job}', wd = tmp_path)

    """Create unstructured grid data."""
    with h5py.File(tmp_path/f'{job}.hdf5', 'r') as f:

        PETSc_version = f'{f.attrs["PETSc_version_major"]}.{f.attrs["PETSc_version_minor"]}.{f.attrs["PETSc_version_subminor"]}'
        if damask.util.version(PETSc_version) < '3.24.1': pytest.xfail()

        x_n = f['/geometry']['x_n']                                                                 # nodes coordinates
        T_e = f['/geometry']['T_e']                                                                 # element connectivity

        assert T_e.shape[1] == nodes[polytope][int(order)-1]                                        # nodes per cell

        nEls = T_e.shape[0]
        npc  = T_e.shape[1]
        nPts = x_n.shape[0]

        pts  = vtkPoints()
        pts.SetData(np_to_vtk(np.ascontiguousarray(x_n)))

        cells = vtkCellArray()
        cells.SetNumberOfCells(nEls)
        nT_e = np.concatenate((np.ones((nEls,1),dtype=np.int64)*npc,T_e),axis=1).ravel()
        cells.SetCells(nEls, np_to_vtkId(nT_e))

    if (polytope == 'TRI'):
        el = vtkLagrangeTriangle()
    elif (polytope == 'TET'):
        el = vtkLagrangeTetra()
    uGrid = vtkUnstructuredGrid()
    uGrid.SetPoints(pts)
    uGrid.SetCells(el.GetCellType(), cells)

    """Assert areas/volumes."""
    sizeFilter = vtkCellSizeFilter()
    sizeFilter.SetInputData(uGrid)
    sizeFilter.Update()

    geom = 'Area' if uGrid.GetCell(0).GetCellDimension() == 2 else 'Volume'
    vals = sizeFilter.GetOutput().GetCellData().GetArray(geom)
    size = sum([vals.GetValue(i) for i in range(vals.GetNumberOfTuples())])
    assert np.isclose(size, sizes[polytope])
