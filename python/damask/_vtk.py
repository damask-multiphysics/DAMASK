# SPDX-License-Identifier: AGPL-3.0-or-later
import contextlib
import logging
import os
import threading
from os import PathLike
from pathlib import Path
from typing import Any, Literal

import numpy as np
import numpy.typing as npt

from ._typehints import StrSequence


# needed for visualization but might not be available everywhere
# https://gitlab.kitware.com/vtk/vtk/-/issues/19687
with contextlib.suppress(ImportError):
    import vtkmodules.vtkRenderingOpenGL2                                                           # noqa

from vtkmodules.vtkCommonCore import (
    vtkVersion,
    vtkPoints,
    vtkStringArray,
    vtkLookupTable,
)

if has_vtkhdf := (np.lib.NumpyVersion(vtkVersion.GetVTKVersion()) >= '9.4.0'):
    from vtkmodules.vtkIOHDF import vtkHDFReader, vtkHDFWriter

from vtkmodules.vtkCommonDataModel import (
    vtkDataSet,
    vtkCellArray,
    vtkImageData,
    vtkRectilinearGrid,
    vtkUnstructuredGrid,
    vtkPointData,
    vtkCellData,
    vtkPolyData,
    VTK_LAGRANGE_TRIANGLE,
    VTK_LAGRANGE_QUADRILATERAL,
    VTK_LAGRANGE_TETRAHEDRON,
    VTK_LAGRANGE_HEXAHEDRON,
)
from vtkmodules.vtkIOLegacy import (
    vtkGenericDataObjectReader,
    vtkDataSetWriter,
)
from vtkmodules.vtkIOXML import (
    vtkXMLWriter,
    vtkXMLImageDataReader,
    vtkXMLImageDataWriter,
    vtkXMLRectilinearGridReader,
    vtkXMLRectilinearGridWriter,
    vtkXMLUnstructuredGridReader,
    vtkXMLUnstructuredGridWriter,
    vtkXMLPolyDataReader,
    vtkXMLPolyDataWriter,
)
from vtkmodules.vtkRenderingCore import (
    vtkDataSetMapper,
    vtkActor,
    vtkRenderer,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
)
from vtkmodules.vtkRenderingAnnotation import (
    vtkScalarBarActor,
)
from vtkmodules.util.numpy_support import (
    numpy_to_vtk,
    numpy_to_vtkIdTypeArray,
    vtk_to_numpy,
)
from vtkmodules.vtkFiltersVerdict import (
    vtkCellSizeFilter
)

from . import Colormap, Table, util
from ._typehints import FloatSequence, IntSequence


logger = logging.getLogger(__name__)

class VTK:
    """
    Spatial visualization (and potentially manipulation).

    High-level interface to VTK.
    """

    def __init__(self,
                 vtk_data: vtkDataSet):
        """
        New spatial visualization.

        Parameters
        ----------
        vtk_data : subclass of vtkDataSet
            Description of geometry and topology, optionally with attached data.
            Valid types are vtkImageData, vtkUnstructuredGrid,
            vtkPolyData, and vtkRectilinearGrid.
        """
        self.vtk_data: vtkDataSet = vtk_data


    def __repr__(self) -> str:
        """
        Return repr(self).

        Give short, human-readable summary.
        """
        info = [self.vtk_data.__vtkname__]                                                          # type: ignore[attr-defined]

        for data in ['Cell Data', 'Point Data']:
            if data == 'Cell Data':  info.append(f'\n# cells: {self.N_cells}')
            if data == 'Point Data': info.append(f'\n# points: {self.N_points}')
            if data in self.labels:
                info += [f'  - {l}' for l in self.labels[data]]

        return util.srepr(info)


    def __eq__(self,
               other: object) -> bool:
        """
        Return self==other.

        Test equality of other.

        Parameters
        ----------
        other : damask.VTK
            VTK to check for equality.

        Returns
        -------
        equal : bool
            Whether both arguments are equal.
        """
        if not isinstance(other, VTK):
            return NotImplemented
        return self.as_ASCII() == other.as_ASCII()



    def copy(self):
        """Return a deep copy."""
        match self.vtk_data:
            case vtkImageData():
                dup = vtkImageData()
            case vtkUnstructuredGrid():
                dup = vtkUnstructuredGrid()
            case vtkPolyData():
                dup = vtkPolyData()
            case vtkRectilinearGrid():
                dup = vtkRectilinearGrid()

        dup.DeepCopy(self.vtk_data)

        return VTK(dup)


    @property
    def comments(self) -> list[str]:
        """Comments in vtkdata"""                                                                   # noqa: D415
        field_data = self.vtk_data.GetFieldData()
        for a in range(field_data.GetNumberOfArrays()):
            if field_data.GetArrayName(a) == 'comments':
                comments = field_data.GetAbstractArray(a)
                return [comments.GetValue(i) for i in range(comments.GetNumberOfValues())]
        return []

    @comments.setter
    def comments(self,
                 comments: StrSequence):
        """
        Set comments.

        Parameters
        ----------
        comments : sequence of str
            Comments to assign to vtkdata.
        """
        s = vtkStringArray()
        s.SetName('comments')
        for c in comments:
            s.InsertNextValue(c)
        self.vtk_data.GetFieldData().AddArray(s)


    @property
    def N_points(self) -> int:
        """Number of points in vtkdata"""                                                           # noqa: D415
        return self.vtk_data.GetNumberOfPoints()


    @property
    def N_cells(self) -> int:
        """Number of cells in vtkdata"""                                                            # noqa: D415
        return self.vtk_data.GetNumberOfCells()


    @property
    def labels(self):
        """Labels of datasets"""                                                                    # noqa: D415
        labels = {}

        cell_data = self.vtk_data.GetCellData()
        if c := [cell_data.GetArrayName(a) for a in range(cell_data.GetNumberOfArrays())]:
            labels['Cell Data'] = c

        point_data = self.vtk_data.GetPointData()
        if p := [point_data.GetArrayName(a) for a in range(point_data.GetNumberOfArrays())]:
            labels['Point Data'] = p

        return labels


    @property
    def measures(self) -> npt.NDArray[np.floating] | None:
        """Areas/volumes of all cells in vtkdata"""                                                 # noqa: D415
        if isinstance(self.vtk_data,vtkUnstructuredGrid):
            size_filter = vtkCellSizeFilter()
            size_filter.SetInputData(self.vtk_data)
            size_filter.Update()
            m = size_filter.GetOutput().GetCellData().GetArray( \
                    'Area' if self.vtk_data.GetCell(0).GetCellDimension() == 2 else 'Volume')
            measures = np.array([m.GetValue(i) for i in range(m.GetNumberOfTuples())])
        else:
            measures = None

        return measures


    @classmethod
    def from_image_data(cls,
                        cells: IntSequence,
                        size: FloatSequence,
                        origin: FloatSequence = np.zeros(3)) -> 'VTK':
        """
        Create VTK of type vtkImageData.

        This is the common type for grid solver results.

        Parameters
        ----------
        cells : sequence of int, len (3)
            Number of cells along each dimension.
        size : sequence of float, len (3)
            Edge length along each dimension.
        origin : sequence of float, len (3), optional
            Coordinates of grid origin.

        Returns
        -------
        new : damask.VTK
            VTK-based geometry without nodal or cell data.

        Examples
        --------
        Create image data with larger spacing along z-direction:

        >>> import damask
        >>> cells = (16,8,4)
        >>> size = (1.0,0.5,.4)
        >>> print(v := damask.VTK.from_image_data(cells=cells,size=size))
        vtkImageData
        <BLANKLINE>
        # cells: 512
        <BLANKLINE>
        # points: 765
        """
        vtk_data = vtkImageData()
        vtk_data.SetDimensions(*(np.array(cells)+1))
        vtk_data.SetOrigin(*(np.array(origin)))
        vtk_data.SetSpacing(*(np.array(size)/np.array(cells)))

        return cls(vtk_data)


    @classmethod
    def from_unstructured_grid(cls,
                               nodes: npt.NDArray[np.floating],
                               connectivity: npt.NDArray[np.integer],
                               cell_type: str,
                               points: npt.NDArray[np.floating] | None = None,
                               extrapolate: Literal[None, 'average', 'weighted'] = None
                              ) -> 'VTK':
        """
        Create VTK of type vtkUnstructuredGrid.

        This is the common type for mesh solver results.
        Elements are of the Lagrange type and all elements
        should have the same type and order.

        Parameters
        ----------
        nodes : numpy.ndarray, shape (:,3)
            Spatial position of the nodes.
        connectivity : numpy.ndarray of np.dtype = np.int64
            Cell connectivity (0-based), first dimension determines #Cells,
            second dimension determines #Nodes/Cell.
        cell_type : str
            Name of the vtkCell subclass. Accepts variation with or without the
            full VTK cell name e.g. "QUAD", "vtk_lagrange_quadrilateral", etc. are
            all valid options.
        points: numpy.ndarray, shape(:,3), optional
            Spatial position of the integration points.
        extrapolate: str, optional
            Extrapolation averaging method.

        Returns
        -------
        new : damask.VTK
            VTK-based geometry without nodal or cell data.

        Examples
        --------
        Create a first-order tetrahedron:

        >>> import damask
        >>> import numpy as np
        >>> nodes = np.array([[0,0,0],[1,0,0],[1,1,0],[0,0,1]])
        >>> connectivity = np.array([[0,1,2,3]])
        >>> print(v := damask.VTK.from_unstructured_grid(nodes=nodes,connectivity=connectivity,
        ...                                              cell_type='TETRAHEDRON'))
        vtkUnstructuredGrid
        <BLANKLINE>
        # cells: 1
        <BLANKLINE>
        # points: 4
        """

        def quad_eq_area_node(pts: npt.NDArray[np.floating]
                              ) -> npt.NDArray[np.floating]:
            """Node that divides a quadrilateral into four equal area quadrilaterals."""
            p_0, p_1, p_2, p_3 = pts[0], pts[1], pts[2], pts[3]
            d = (p_2[0]-p_0[0])*(p_3[1]-p_1[1])+(p_1[0]-p_3[0])*(p_2[1]-p_0[1])
            x = 0.5*(p_0[0]+p_1[0]+p_2[0]+p_3[0]) - \
                ((p_1[0]*p_2[0]-p_0[0]*p_1[0])*p_3[1]+(p_0[0]*p_1[0]-p_0[0]*p_3[0])*p_2[1] + \
                 (p_0[0]-p_2[0])*p_3[0]*p_1[1]+(p_2[0]*p_3[0]-p_1[0]*p_2[0])*p_0[1]) / d
            y = 0.5*(p_0[1]+p_1[1]+p_2[1]+p_3[1]) - \
                (((p_1[0]-p_0[0])*p_2[1]+(p_2[0]-p_1[0])*p_0[1])*p_3[1] + \
                 ((p_0[0]-p_3[0])*p_2[1]+(p_3[0]-p_2[0])*p_0[1])*p_1[1]) / d

            return np.array([x, y, 0.0])

        def split_elements(nodes: npt.NDArray[np.floating],
                           points: npt.NDArray[np.floating],
                           connectivity: npt.NDArray[np.integer],
                           cell_type: str,
                           ) -> tuple[npt.NDArray[np.floating], npt.NDArray[np.integer]]:
            """Split n-th order elements down to first order elements."""
            n_cells = connectivity.shape[0]
            n_nodes = nodes.shape[0]
            n_ips   = points.shape[0] // n_cells    # per cell

            nodes_  = T_e = None
            idx     = np.empty(0, dtype=np.int64)

            if cell_type == 'TRIA':
                match n_ips:
                    case  1:                                                                        # quad = 'o1f'
                        pass
                    case  3:                                                                        # quad = 'o2f'
                        idx = np.array([[0,3,4],[0,4,5],                                            # IP0
                                        [5,4,2],                                                    # IP1
                                        [3,1,4]])                                                   # IP2
                    case  6:                                                                        # quad = 'o3f'
                        idx = np.array([[5,6,9],                                                    # IP0
                                        [0,3,9],[0,9,8],                                            # IP1
                                        [3,4,5],[3,5,9],                                            # IP2
                                        [2,7,6],                                                    # IP3
                                        [6,8,9],[6,7,8],                                            # IP4
                                        [1,5,4]])                                                   # IP5
                    case  7:                                                                        # quad = 'o4r'
                        idx = np.array([[3,13,12],[12,13,14],[11,12,14],                            # IP0
                                        [0,3,12],[0,12,11],                                         # IP1
                                        [6,7,13],[7,14,13],[7,8,14],                                # IP2
                                        [2,9,8],[8,9,14],                                           # IP3
                                        [3,4,13],[4,5,13],                                          # IP4
                                        [1,6,5],[5,6,13],                                           # IP5
                                        [10,11,14],[9,10,14]])                                      # IP6
                    case 12:                                                                        # quad = 'o4f'
                        idx = np.array([[0,3,12],[0,12,11],                                         # IP0
                                        [3,13,12],[11,12,14],                                       # IP1
                                        [2,9,8],                                                    # IP2,
                                        [7,14,12],                                                  # IP3
                                        [1,6,5],                                                    # IP4
                                        [7,12,13],                                                  # IP5
                                        [8,10,14],[8,9,10],                                         # IP6
                                        [10,11,14],                                                 # IP7
                                        [4,5,6],[4,6,13],                                           # IP8
                                        [3,4,13],                                                   # IP9
                                        [7,8,14],                                                   # IP10
                                        [6,7,13]])                                                  # IP11
                    case 15:                                                                        # quad = 'o5r'
                        idx = np.array([[0,3,15],[0,15,14],                                         # IP0
                                        [15,18,19],[15,19,20],                                      # IP1
                                        [8,9,19],                                                   # IP2
                                        [2,11,10],                                                  # IP3
                                        [9,20,19],[9,17,20],                                        # IP4
                                        [4,5,16],[4,16,18],                                         # IP5
                                        [1,7,6],                                                    # IP6
                                        [8,18,16],[8,19,18],                                        # IP7
                                        [13,20,17],[12,13,17],                                      # IP8
                                        [10,12,17],[10,11,12],                                      # IP9
                                        [14,15,20],[13,14,20],                                      # IP10
                                        [5,6,7],[5,7,16],                                           # IP11
                                        [3,4,18],[3,18,15],                                         # IP12
                                        [9,10,17],                                                  # IP13
                                        [7,8,16]])                                                  # IP14
                    case 16:                                                                        # quad = 'o5f'
                        idx = np.array([[16,19,18],[18,19,20],[17,20,19],                           # IP0
                                        [8,19,16],[8,9,19],[9,17,19],                               # IP1
                                        [3,15,14],[3,18,15],[15,18,20],[14,15,20],                  # IP2
                                        [0,3,14],                                                   # IP3
                                        [4,5,16],[4,16,18],                                         # IP4
                                        [10,12,17],                                                 # IP5
                                        [2,11,10],                                                  # IP6
                                        [13,20,17],[12,13,17],                                      # IP7
                                        [5,7,16],                                                   # IP8
                                        [1,7,6],                                                    # IP9
                                        [10,11,12],                                                 # IP10
                                        [13,14,20],                                                 # IP11
                                        [5,6,7],                                                    # IP12
                                        [3,4,18],                                                   # IP13
                                        [9,10,17],                                                  # IP14
                                        [7,8,16]])                                                  # IP15
                    case  _:
                        raise ValueError('visualization not implemented for triangles with p_s > 5')

                T_e = np.vstack([connectivity[k,idx] for k in range(n_cells)]) if len(idx) else \
                      connectivity

            elif cell_type == 'QUAD':
                match n_ips:
                    case  1:                                                                        # quad = 'o1r'
                        pass
                    case  4:                                                                        # quad = 'o1f'
                        ns = lambda k,l : np.array(np.sum(nodes[connectivity[k,np.array(l)]],axis=0)/len(l))
                        nodes_ = np.vstack(([(ns(k,[0,1]), ns(k,[1,2]), ns(k,[2,3]), ns(k,[3,0]),
                                              quad_eq_area_node(nodes[connectivity[k,[0,1,2,3]]]))
                                            for k in range(n_cells)]))
                        n5p = lambda k,i: n_nodes+(9-4)*k+i
                        T_e = np.array([(connectivity[k,0], n5p(k,0), n5p(k,4), n5p(k,3),           # IP0
                                         n5p(k,0), connectivity[k,1], n5p(k,1), n5p(k,4),           # IP2
                                         n5p(k,3), n5p(k,4), n5p(k,2), connectivity[k,3],           # IP1
                                         n5p(k,4), n5p(k,1), connectivity[k,2], n5p(k,2))           # IP3
                                        for k in range(n_cells)]).reshape(-1,4)
                    case  9:                                                                        # quad = 'o2f'
                        nd = lambda k,i,j,w : (nodes[connectivity[k,i]] + w \
                                              * (nodes[connectivity[k,j]]-nodes[connectivity[k,i]]))
                        nodes_ = np.vstack([(n_09 := nd(k,0,1,1/3), n_10 := nd(k,0,1,2/3),
                                             nd(k,1,2,1/3), nd(k,1,2,2/3),
                                             n_13 := nd(k,2,3,1/3), n_14 := nd(k,2,3,2/3),
                                             nd(k,3,0,1/3), nd(k,3,0,2/3),
                                             n_09+1/3*(n_14-n_09), n_10+1/3*(n_13-n_10),
                                             n_10+2/3*(n_13-n_10), n_09+2/3*(n_14-n_09))
                                            for k in range(n_cells)])
                        n12p = lambda k,l : n_nodes+(16-4)*k+np.array(l)
                        T_e = np.vstack([np.concatenate(
                                          ([connectivity[k,0]], n12p(k,[0,8,7]),                    # IP 0
                                           n12p(k,[0,1,9,8]),                                       # IP 3
                                           n12p(k,[1]), [connectivity[k,1]], n12p(k,[2,9]),         # IP 6
                                           n12p(k,[7,8,11,6]),                                      # IP 1
                                           n12p(k,[8,9,10,11]),                                     # IP 4
                                           n12p(k,[9,2,3,10]),                                      # IP 7
                                           n12p(k,[6,11,5]), [connectivity[k,3]],                   # IP 2
                                           n12p(k,[11,10,4,5]),                                     # IP 5
                                           n12p(k,[10,3]),[connectivity[k,2]],n12p(k,[4]))          # IP 8
                                         ).reshape(-1,4) for k in range(n_cells)])
                    case 16:                                                                        # quad = 'o3f'
                        nd = lambda k,i,j,w : (nodes[connectivity[k,i]] + w \
                                              * (nodes[connectivity[k,j]]-nodes[connectivity[k,i]]))
                        nodes_ = np.vstack([(n_17 := nd(k,0,1,1/4), n_18 := nd(k,0,1,2/4),
                                             n_19 := nd(k,0,1,3/4),
                                             nd(k,1,2,1/4), nd(k,1,2,2/4), nd(k,1,2,3/4),
                                             n_23 := nd(k,2,3,1/4), n_24 := nd(k,2,3,2/4),
                                             n_25 := nd(k,2,3,3/4),
                                             nd(k,3,0,1/4), nd(k,3,0,2/4), nd(k,3,0,3/4),
                                             n_17+1/4*(n_25-n_17), n_19+1/4*(n_23-n_19),
                                             n_19+3/4*(n_23-n_19), n_17+3/4*(n_25-n_17),
                                             n_18+1/4*(n_24-n_18), n_19+2/4*(n_23-n_19),
                                             n_18+3/4*(n_24-n_18), n_17+2/4*(n_25-n_17),
                                             n_18+2/4*(n_24-n_18))
                                            for k in range(n_cells)])
                        n21p = lambda k, l : n_nodes+(25-4)*k+np.array(l)
                        T_e = np.vstack([np.concatenate(
                                          ([connectivity[k,0]], n21p(k,[0,12,11]),                  # IP 0
                                           n21p(k,[0,1,16,12]),                                     # IP 4
                                           n21p(k,[1,2,13,16]),                                     # IP 8
                                           n21p(k,[2]), [connectivity[k,1]], n21p(k,[3,13]),        # IP 12
                                           n21p(k,[11,12,19,10]),                                   # IP 1
                                           n21p(k,[12,16,20,19]),                                   # IP 5
                                           n21p(k,[16,13,17,20]),                                   # IP 9
                                           n21p(k,[13,3,4,17]),                                     # IP 13
                                           n21p(k,[10,19,15,9]),                                    # IP 2
                                           n21p(k,[19,20,18,15]),                                   # IP 6
                                           n21p(k,[20,17,14,18]),                                   # IP 10
                                           n21p(k,[17,4,5,14]),                                     # IP 14
                                           n21p(k,[9,15,8]), [connectivity[k,3]],                   # IP 3
                                           n21p(k,[15,18,7,8]),                                     # IP 7
                                           n21p(k,[18,14,6,7]),                                     # IP 11
                                           n21p(k,[14,5]), [connectivity[k,2]], n21p(k,[6]))        # IP 15
                                         ).reshape(-1,4) for k in range(n_cells)])
                    case  _:
                        raise ValueError('visualization not implemented for quadrilaterals with p_s > 3')

            elif cell_type == 'TETRA':
                match n_ips:
                    case 1:                                                                         # quad = 'o1f'
                        pass
                    case 4:                                                                         # quad = 'o2f'
                        idx = np.array([[6,5,2,9],[0,5,6,9],                                        # IP0
                                        [1,5,4,8],[0,4,5,8],                                        # IP1
                                        [0,8,5,9],[0,7,8,9],                                        # IP2
                                        [7,8,9,3]])                                                 # IP3
                    case _:
                        raise ValueError('visualization not implemented for tetrahedra with p_s > 2')

                T_e = np.vstack([connectivity[k,idx] for k in range(n_cells)]) if len(idx) else \
                      connectivity

            elif cell_type == 'HEXA':
                match n_ips:
                    case  1:                                                                        # quad= 'o1r'
                        pass
                    case  8:                                                                        # quad= 'o1f'
                        ns = lambda k,l : np.sum(nodes[connectivity[k,l]],axis=0)/len(l)
                        nodes_ = np.vstack([(ns(k,[0,1]), ns(k,[1,2]), ns(k,[2,3]), ns(k,[3,0]),
                                             ns(k,[4,5]), ns(k,[5,6]), ns(k,[6,7]), ns(k,[7,4]),
                                             ns(k,[0,4]), ns(k,[1,5]), ns(k,[2,6]), ns(k,[3,7]),
                                             ns(k,[0,1,4,5]), ns(k,[1,2,5,6]), ns(k,[2,3,6,7]),
                                             ns(k,[0,3,4,7]), ns(k,[0,1,2,3]), ns(k,[4,5,6,7]),
                                             ns(k,np.arange(8)))
                                            for k in range(n_cells)])
                        n19p = lambda k,l : n_nodes+(27-8)*k+np.array(l)
                        T_e = np.vstack([np.concatenate(
                                          ([connectivity[k,0]], n19p(k,[0,16,3,8,12,18,15]),        # IP 0
                                           n19p(k,[0]), [connectivity[k,1]], n19p(k,[1,16,12,9,13,18]),
                                           n19p(k,[16,1]), [connectivity[k,2]], n19p(k,[2,18,13,10,14]),
                                           n19p(k,[3,16,2]), [connectivity[k,3]], n19p(k,[15,18,14,11]),
                                           n19p(k,[8,12,18,15]), [connectivity[k,4]], n19p(k,[4,17,7]),
                                           n19p(k,[12,9,13,18,4]), [connectivity[k,5]], n19p(k,[5,17]),
                                           n19p(k,[18,13,10,14,17,5]), [connectivity[k,6]], n19p(k,[6]),
                                           n19p(k,[15,18,14,11,7,17,6]), [connectivity[k,7]])       # IP 3
                                          ).reshape(-1,8) for k in range(n_cells)])
                    case 27:                                                                        # quad= 'o2f'
                        nd = lambda k,i,j,w : (nodes[connectivity[k,i]] + w \
                                              * (nodes[connectivity[k,j]]-nodes[connectivity[k,i]]))
                        nt = lambda i,j,w : i + w * (j - i)
                        nodes_ = np.vstack([np.vstack(
                                            (n_27 := nd(k,0,1,1/3), n_28 := nd(k,0,1,2/3),
                                             nd(k,1,2,1/3), nd(k,1,2,2/3),
                                             n_31 := nd(k,2,3,1/3), n_32 := nd(k,2,3,2/3),
                                             nd(k,3,0,1/3), nd(k,3,0,2/3),
                                             nt(n_27,n_32,1/3), nt(n_28,n_31,1/3),
                                             nt(n_28,n_31,2/3), nt(n_27,n_32,2/3),
                                             n_39:= nd(k,0,4,1/3), n_40 := nd(k,1,5,1/3),
                                             n_41:= nd(k,2,6,1/3), n_42 := nd(k,3,7,1/3),
                                             n_43 := nt(n_39,n_40,1/3), n_44 := nt(n_39,n_40,2/3),
                                             nt(n_40,n_41,1/3), nt(n_40,n_41,2/3),
                                             n_47 := nt(n_41,n_42,1/3), n_48 := nt(n_41,n_42,2/3),
                                             nt(n_42,n_39,1/3), nt(n_42,n_39,2/3),
                                             nt(n_43,n_48,1/3), nt(n_44,n_47,1/3),
                                             nt(n_44,n_47,2/3), nt(n_43,n_48,2/3),
                                             n_55 := nd(k,0,4,2/3), n_56 := nd(k,1,5,2/3),
                                             n_57 := nd(k,2,6,2/3), n_58 := nd(k,3,7,2/3),
                                             n_59 := nt(n_55,n_56,1/3), n_60 := nt(n_55,n_56,2/3),
                                             nt(n_56,n_57,1/3), nt(n_56,n_57,2/3),
                                             n_63 := nt(n_57,n_58,1/3), n_64 := nt(n_57,n_58,2/3),
                                             nt(n_58,n_55,1/3), nt(n_58,n_55,2/3),
                                             nt(n_59,n_64,1/3), nt(n_60,n_63,1/3),
                                             nt(n_60,n_63,2/3), nt(n_59,n_64,2/3),
                                             n_71 := nd(k,4,5,1/3), n_72 := nd(k,4,5,2/3),
                                             nd(k,5,6,1/3), nd(k,5,6,2/3),
                                             n_75 := nd(k,6,7,1/3), n_76 := nd(k,6,7,2/3),
                                             nd(k,7,4,1/3), nd(k,7,4,2/3),
                                             nt(n_71,n_76,1/3), nt(n_72,n_75,1/3),
                                             nt(n_72,n_75,2/3), nt(n_71,n_76,2/3)))
                                            for k in range(n_cells)])
                        n56p = lambda k, l : n_nodes+(64-8)*k+np.array(l)
                        T_e = np.vstack([np.concatenate(
                                         ([connectivity[k,0]], n56p(k,[0,8,7,12,16,24,23]),         # IP 0
                                          n56p(k,[0,1,9,8,16,17,25,24]),                            # IP 9
                                          n56p(k,[1]), [connectivity[k,1]], n56p(k,[2,9,17,+13,+18,25]),
                                          n56p(k,[7,8,11,6,23,24,27,22]),
                                          n56p(k,[8,9,10,11,24,25,26,27]),
                                          n56p(k,[9,2,3,10,25,18,19,26]),
                                          n56p(k,[6,11,5]), [connectivity[k,3]], n56p(k,[22,27,21,15]),
                                          n56p(k,[11,10,4,5,27,26,20,21]),
                                          n56p(k,[10,3]), [connectivity[k,2]], n56p(k,[4,26,19,14,20]),
                                          n56p(k,[12,16,24,23,28,32,40,39]),                        # IP 1
                                          n56p(k,[16,17,25,24,32,33,41,40]),                        # IP 10
                                          n56p(k,[17,13,18,25,33,29,34,41]),                        # IP 19
                                          n56p(k,[23,24,27,22,39,40,43,38]),
                                          n56p(k,[24,25,26,27,40,41,42,43]),
                                          n56p(k,[25,18,19,26,41,34,35,42]),
                                          n56p(k,[22,27,21,15,38,43,37,31]),
                                          n56p(k,[27,26,20,21,43,42,36,37]),
                                          n56p(k,[26,19,14,20,42,35,30,36]),
                                          n56p(k,[28,32,40,39]), [connectivity[k,4]], n56p(k,[44,52,51]),
                                          n56p(k,[32,33,41,40,44,45,53,52]),
                                          n56p(k,[33,29,34,41,45]), [connectivity[k,5]], n56p(k,[46,53]),
                                          n56p(k,[39,40,43,38,51,52,55,50]),
                                          n56p(k,[40,41,42,43,52,53,54,55]),
                                          n56p(k,[41,34,35,42,53,46,47,54]),
                                          n56p(k,[38,43,37,31,50,55,49]), [connectivity[k,7]],      # IP 8
                                          n56p(k,[43,42,36,37,55,54,48,49]),                        # IP 17
                                          n56p(k,[42,35,30,36,54,47]), [connectivity[k,6]] ,n56p(k,[48]))
                                         ).reshape(-1,8) for k in range(n_cells)])
                    case  _:
                        raise ValueError('visualization not implemented for hexahedra with p_s > 2')

            nodes_ = nodes if nodes_ is None else np.vstack((nodes, nodes_))
            T_e    = connectivity if T_e is None else T_e

            return nodes_, T_e

        cell_types = {'TRIA':VTK_LAGRANGE_TRIANGLE,    'QUAD':VTK_LAGRANGE_QUADRILATERAL,
                      'TETR':VTK_LAGRANGE_TETRAHEDRON, 'HEXA':VTK_LAGRANGE_HEXAHEDRON}
        vtk_cell_type = cell_type.split('_')[-1].upper()[:4]
        if vtk_cell_type not in cell_types.keys():
            raise KeyError(f'unknown cell type {cell_type}')

        if points is not None and extrapolate is None:
            nodes, connectivity = split_elements(nodes, points, connectivity, vtk_cell_type)

        vtk_nodes = vtkPoints()
        vtk_nodes.SetData(numpy_to_vtk(np.ascontiguousarray(nodes)))
        vtk_cells = vtkCellArray()
        vtk_cells.SetData(connectivity.shape[1],numpy_to_vtkIdTypeArray(connectivity.astype(np.int64).ravel()))
        vtk_data = vtkUnstructuredGrid()
        vtk_data.SetPoints(vtk_nodes)
        vtk_data.SetCells(cell_types[vtk_cell_type],vtk_cells)

        return cls(vtk_data)


    @classmethod
    def from_poly_data(cls,
                       points: npt.NDArray[np.floating]) -> 'VTK':
        """
        Create VTK of type polyData.

        This is the common type for point-wise data.

        Parameters
        ----------
        points : numpy.ndarray, shape (:,3)
            Spatial position of the points.

        Returns
        -------
        new : damask.VTK
            VTK-based geometry without nodal or cell data.
        """
        N = points.shape[0]
        vtk_points = vtkPoints()
        vtk_points.SetData(numpy_to_vtk(np.ascontiguousarray(points)))

        vtk_cells = vtkCellArray()
        vtk_cells.SetNumberOfCells(N)
        vtk_cells.SetCells(N,numpy_to_vtkIdTypeArray(np.stack((np.ones  (N,dtype=np.int64),
                                                            np.arange(N,dtype=np.int64)),axis=1).ravel(),deep=True))

        vtk_data = vtkPolyData()
        vtk_data.SetPoints(vtk_points)
        vtk_data.SetVerts(vtk_cells)

        return cls(vtk_data)


    @classmethod
    def from_rectilinear_grid(cls,coordinates: FloatSequence) -> 'VTK':
        """
        Create VTK of type vtkRectilinearGrid.

        Parameters
        ----------
        coordinates : sequence of sequences of floats, len (3)
            Grid coordinates along x, y, and z directions.

        Returns
        -------
        new : damask.VTK
            VTK-based geometry without nodal or cell data.
        """
        vtk_data = vtkRectilinearGrid()
        vtk_data.SetDimensions(*map(len,coordinates))
        coord = [numpy_to_vtk(np.array(coordinates[i]),deep=True) for i in [0,1,2]]
        [coord[i].SetName(n) for i,n in enumerate(['x','y','z'])]
        vtk_data.SetXCoordinates(coord[0])
        vtk_data.SetYCoordinates(coord[1])
        vtk_data.SetZCoordinates(coord[2])

        return cls(vtk_data)


    @classmethod
    def load(cls,
             fname: str | PathLike,
             dataset_type: Literal[None, 'ImageData', 'UnstructuredGrid', 'PolyData', 'RectilinearGrid'] = None) -> 'VTK':
        """
        Load from VTK file.

        Parameters
        ----------
        fname : str or pathlib.Path
            Filename to read.
            Valid extensions are .vti, .vtu, .vtp, .vtr, vtkhdf, and .vtk.
        dataset_type : {'ImageData', 'UnstructuredGrid', 'PolyData', 'RectilinearGrid'}, optional
            Name of the vtkDataSet subclass when opening a .vtk file.

        Returns
        -------
        loaded : damask.VTK
            VTK-based geometry from file.

        Notes
        -----
        Loading VTKHDF files requires VTK 9.4 or later for PolyData, UnstructuredGrid,
        and ImageData and VTK 9.7 or later for RectilinearGrid.
        """
        if not Path(fname).expanduser().is_file():                                                  # vtk has a strange error handling
            raise FileNotFoundError(f'file "{fname}" not found')

        if (ext := Path(fname).suffix) == '.vtk' or dataset_type is not None:
            reader_legacy = vtkGenericDataObjectReader()
            reader_legacy.SetFileName(str(Path(fname).expanduser()))
            reader_legacy.Update()
            if dataset_type is None:
                raise TypeError('missing dataset type for legacy VTK file')
            dtl = dataset_type.lower()
            vtk_data = (
                reader_legacy.GetStructuredPointsOutput() if dtl.endswith(('imagedata', 'image_data')) else
                reader_legacy.GetUnstructuredGridOutput() if dtl.endswith(('unstructuredgrid', 'unstructured_grid')) else
                reader_legacy.GetPolyDataOutput() if dtl.endswith(('polydata', 'poly_data')) else
                reader_legacy.GetRectilinearGridOutput() if dtl.endswith(('rectilineargrid', 'rectilinear_grid')) else
                None
            )
            if vtk_data is None:
                raise TypeError(f'unsupported VTK dataset type "{dataset_type}"')
        else:
            reader = (
                vtkXMLImageDataReader() if ext == '.vti' else
                vtkXMLUnstructuredGridReader() if ext == '.vtu' else
                vtkXMLPolyDataReader() if ext == '.vtp' else
                vtkXMLRectilinearGridReader() if ext == '.vtr' else
                vtkHDFReader() if (ext == '.vtkhdf' and has_vtkhdf) else
                None
            )
            if reader is None:
                raise TypeError(f'unsupported VTK file extension "{ext}"')
            reader.SetFileName(str(Path(fname).expanduser()))
            reader.Update()
            vtk_data = reader.GetOutputAsDataSet()

        return cls(vtk_data)


    def as_ASCII(self) -> str:
        """ASCII representation of the VTK data."""
        writer = vtkDataSetWriter()
        writer.SetHeader(f'# {util.execution_stamp("VTK")}')
        writer.WriteToOutputStringOn()
        writer.SetInputData(self.vtk_data)
        writer.Write()
        return writer.GetOutputString()


    @staticmethod
    def _write(writer):
        """Wrapper for parallel writing."""
        writer.Write()

    def save(self,
             fname: str | PathLike,
             parallel: bool = True,
             compress: bool = True):
        """
        Save as VTK file.

        Parameters
        ----------
        fname : str or pathlib.Path
            Filename to write.
        parallel : bool, optional
            Write data in parallel background process. Defaults to True.
        compress : bool, optional
            Compress with zlib algorithm. Defaults to True.
        """
        writer: vtkXMLWriter | None = (
            vtkXMLImageDataWriter() if isinstance(self.vtk_data, vtkImageData) else
            vtkXMLUnstructuredGridWriter() if isinstance(self.vtk_data, vtkUnstructuredGrid) else
            vtkXMLPolyDataWriter() if isinstance(self.vtk_data, vtkPolyData) else
            vtkXMLRectilinearGridWriter() if isinstance(self.vtk_data, vtkRectilinearGrid) else
            None
        )
        assert writer is not None

        default_ext = '.'+writer.GetDefaultFileExtension()
        ext = Path(fname).suffix
        writer.SetFileName(str(Path(fname).expanduser())+(default_ext if default_ext != ext else ''))

        if compress:
            writer.SetCompressorTypeToZLib()
        else:
            writer.SetCompressorTypeToNone()
        writer.SetByteOrderToLittleEndian()
        writer.SetDataModeToBinary()
        writer.SetInputData(self.vtk_data)

        if parallel:
            try:
                writer_th = threading.Thread(target=self._write,args=(writer,))
                writer_th.start()
            except TypeError:
                writer.Write()
        else:
            writer.Write()


    def save_VTKHDF(self,
                    fname: str | PathLike):
        """
        Save as VTKHDF file.

        Parameters
        ----------
        fname : str or pathlib.Path
            Filename to write.

        Notes
        -----
        Saving as VTKHDF file requires VTK 9.4 or later for PolyData,
        and UnstructuredGrid and VTK 9.7 or later for RectilinearGrid and
        ImageData.
        """
        if not has_vtkhdf:
            raise NotImplementedError('save as VTKHDF requires VTK 9.4 or later')

        writer = vtkHDFWriter()
        ext = Path(fname).suffix
        writer.SetFileName(str(Path(fname).expanduser())+('.vtkhdf' if '.vtkhdf' != ext else ''))
        writer.SetInputData(self.vtk_data)
        writer.Write()


    def set(self,
            label: str,
            data: npt.NDArray[Any] | np.ma.MaskedArray,
            component_names: StrSequence | None = None,
            info: str | None = None) -> 'VTK':
        """
        Add new or replace existing point or cell data.

        Parameters
        ----------
        label : str, optional
            Label of data array.
        data : numpy.ndarray or numpy.ma.MaskedArray
            Data to add or replace. First array dimension needs to match either
            number of cells or number of points.
        component_names: sequence of str, optional
            Name of each data component. Flattened component names are assigned to
            flattened shape of one data row, assuming row-major order for both.
        info : str, optional
            Human-readable information about the data.

        Returns
        -------
        updated : damask.VTK
            Updated VTK-based geometry.

        Notes
        -----
        If the number of cells equals the number of points, the data is added to both.


        Examples
        --------
        Add a constant tensor field called 'F' to image data:

        >>> import numpy as np
        >>> import damask
        >>> cells = (16,8,16)
        >>> size = (1.0,0.5,1.0)
        >>> v = damask.VTK.from_image_data(cells=cells,size=size)
        >>> print(v := v.set(label='F',data=np.broadcast_to(np.eye(3),(np.prod(cells),3,3)),
        ...                  component_names = [['11','12','13'],
        ...                                     ['21','22','23'],
        ...                                     ['31','32','33']]))
        vtkImageData
        <BLANKLINE>
        # cells: 2048
          - F
        <BLANKLINE>
        # points: 2601
        """
        dup = self.copy()

        N_p,N_c = dup.vtk_data.GetNumberOfPoints(),dup.vtk_data.GetNumberOfCells()
        if (N_data := data.shape[0]) not in [N_p,N_c]:
            raise ValueError(f'data count mismatch ({N_data} ≠ {N_p} & {N_c})')

        data_ = (data.filled() if isinstance(data,np.ma.MaskedArray) else data) \
                .reshape(N_data,-1) \
                .astype(np.single if data.dtype in [np.double,np.longdouble] else data.dtype)

        if data.dtype.type is np.str_:
            d = vtkStringArray()
            for s in np.squeeze(data_):
                d.InsertNextValue(s)
        else:
            d = numpy_to_vtk(data_,deep=True)

        d.SetName(label)
        if component_names is not None:
            for c,n in zip(range(data_.shape[1]),np.asarray(component_names).flatten(),strict=True):
                d.SetComponentName(c,str(n))

        if N_data == N_p:
            dup.vtk_data.GetPointData().AddArray(d)
        if N_data == N_c:
            dup.vtk_data.GetCellData().AddArray(d)

        if info is not None: dup.comments += [f'{label}: {info}']

        return dup


    def set_from_table(self,
                       table: 'Table',
                       labels: StrSequence | None = None,
                       component_names: dict[str,StrSequence] | None = None,
                       info: str | None = None) -> 'VTK':
        """
        Add new or replace existing point or cell data from damask.Table.

        Parameters
        ----------
        table : damask.Table
            Data to add or replace. Each table label is individually considered.
            Number of rows needs to match either number of cells or number of points.
        labels : sequence of str, optional
            Labels of data to consider. By default, all data is added.
        component_names: dictionary of sequence of str, optional
            Name of each data component per table entry. Flattened component names are
            assigned to flattened data shape of corresponding table entry, assuming
            row-major order for both.
        info : str, optional
            Human-readable information about the data.

        Returns
        -------
        updated : damask.VTK
            Updated VTK-based geometry.

        Notes
        -----
        If the number of cells equals the number of points, the data is added to both.
        """
        dup = self.copy()
        for label in (labels if labels is not None else table.labels):
            dup = dup.set(label,table.get(label),info=info,
                          component_names=(component_names.get(label,None) if component_names else None))

        return dup


    def get(self,
            label: str) -> npt.NDArray[Any]:
        """
        Get either cell or point data.

        Cell data takes precedence over point data, i.e. this
        function assumes that labels are unique among cell and
        point data.

        Parameters
        ----------
        label : str
            Data label.

        Returns
        -------
        data : numpy.ndarray
            Data stored under the given label.
        """
        cell_data: vtkCellData = self.vtk_data.GetCellData()
        point_data: vtkPointData = self.vtk_data.GetPointData()

        if label in [cell_data.GetArrayName(i) for i in range(cell_data.GetNumberOfArrays())]:
            vtk_container: vtkPointData | vtkCellData = cell_data
        elif label in [point_data.GetArrayName(i) for i in range(point_data.GetNumberOfArrays())]:
            vtk_container = point_data
        else:
            raise KeyError(f'array "{label}" not found in cell or point data')

        try:
            return vtk_to_numpy(vtk_container.GetArray(label))
        except AttributeError:                                                                       # assume string array
            vtk_array = vtk_container.GetAbstractArray(label)
            return np.array([vtk_array.GetValue(i) for i in range(vtk_array.GetNumberOfValues())]).astype(str)


    def delete(self,
               label: str) -> 'VTK':
        """
        Delete either cell or point data.

        Cell data takes precedence over point data, i.e. this
        function assumes that labels are unique among cell and
        point data.

        Parameters
        ----------
        label : str
            Data label.

        Returns
        -------
        updated : damask.VTK
            Updated VTK-based geometry.
        """
        dup = self.copy()

        cell_data = dup.vtk_data.GetCellData()
        if label in [cell_data.GetArrayName(a) for a in range(cell_data.GetNumberOfArrays())]:
            dup.vtk_data.GetCellData().RemoveArray(label)
            return dup

        point_data = self.vtk_data.GetPointData()
        if label in [point_data.GetArrayName(a) for a in range(point_data.GetNumberOfArrays())]:
            dup.vtk_data.GetPointData().RemoveArray(label)
            return dup

        raise KeyError(f'array "{label}" not found')


    def show(self,
             label: str | None = None,
             colormap: str | Colormap = 'cividis'):
        """
        Render.

        Parameters
        ----------
        label : str, optional
            Label of the dataset to show.
        colormap : damask.Colormap or str, optional
            Colormap for visualization of dataset. Defaults to 'cividis'.

        Notes
        -----
        The first component is shown when visualizing vector datasets
        (this includes tensor datasets as they are flattened).
        """
        # See http://compilatrix.com/article/vtk-1 for possible improvements.
        if os.name == 'posix' and 'DISPLAY' not in os.environ:
            logger.warning('Found no rendering device')
            return
        try:
            import wx
            _ = wx.App(False)                                                                       # noqa
            width, height = wx.GetDisplaySize()
        except ImportError:
            try:
                import tkinter
                tk = tkinter.Tk()
                width  = tk.winfo_screenwidth()
                height = tk.winfo_screenheight()
                tk.destroy()
            except Exception:
                width  = 1024
                height =  768

        lut = vtkLookupTable()
        colormap_ = Colormap.from_predefined(colormap) if isinstance(colormap,str) else \
                    colormap
        lut.SetNumberOfTableValues(len(colormap_.colors))
        for i,c in enumerate(colormap_.colors):
            lut.SetTableValue(i,c if len(c)==4 else np.append(c,1.0))                               # type: ignore[arg-type]
        lut.Build()
        if label is not None:
            self.vtk_data.GetCellData().SetActiveScalars(label)
        mapper = vtkDataSetMapper()
        mapper.SetInputData(self.vtk_data)
        mapper.SetLookupTable(lut)
        mapper.SetScalarRange(self.vtk_data.GetScalarRange())

        actor = vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(230/255,150/255,68/255)

        ren = vtkRenderer()
        ren.AddActor(actor)
        if label is None:
            ren.SetBackground(67/255,128/255,208/255)
        else:
            colormap_vtk = vtkScalarBarActor()
            colormap_vtk.SetLookupTable(lut)
            colormap_vtk.SetTitle(label)
            colormap_vtk.SetMaximumWidthInPixels(width//100)
            ren.AddActor2D(colormap_vtk)
            ren.SetBackground(0.3,0.3,0.3)

        window = vtkRenderWindow()
        window.AddRenderer(ren)
        window.SetSize(width,height)
        window.SetWindowName(util.execution_stamp('VTK','show'))

        iren = vtkRenderWindowInteractor()
        iren.SetRenderWindow(window)
        window.Render()
        iren.Start()
