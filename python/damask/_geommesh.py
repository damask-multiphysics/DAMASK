# SPDX-License-Identifier: AGPL-3.0-or-later
import itertools
import uuid
from pathlib import Path
from typing import Union

import gmsh
import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy

from . import grid as grid_, util


# (dim, nodes_per_elem) → gmsh element type (linear elements only)
ELEM_TYPES = {
    (2, 3): 2,   # triangle
    (2, 4): 3,   # quadrilateral
    (3, 4): 4,   # tetrahedron
    (3, 8): 5,   # hexahedron
}

# (dim, nodes_per_elem) → gmsh element type for boundary elements
BOUNDARY_TYPES = {
    (0, 1): 15,  # point
    (1, 2): 1,   # line
    (2, 3): 2,   # triangle
    (2, 4): 3,   # quadrilateral
}

# gmsh element type → number of nodes (for reading back from files)
TYPE_TO_NODES = {v: k[1] for k, v in (ELEM_TYPES | BOUNDARY_TYPES).items()}


class GeomMesh:
    """
    Geometry definition for mesh solvers.

    Create and manipulate geometry definitions for storage as gmsh files
    ('.msh' extension). A mesh consists of node coordinates and the
    connectivity for linear elements and contains the material ID
    (indexing an entry in the material configuration). Named boundary
    physical groups store boundary information for applying boundary
    conditions.
    """

    def __init__(self,
                 coordinates: np.ndarray,
                 connectivity: np.ndarray,
                 material: np.ndarray) -> None:
        """
        Create a mesh from node coordinates, element connectivity and material IDs.

        Parameters
        ----------
        coordinates : numpy.ndarray of float, shape (N, dim)
            Node coordinates in meter. dim is 2 or 3.
        connectivity : numpy.ndarray of int, shape (M, nodes_per_elem)
            Element connectivity using 0-based node indices.
        material : numpy.ndarray of int, shape (M,)
            Material indices per element, 0-based.
        """
        if connectivity.shape[0] != material.shape[0]:
            raise ValueError('mismatch between number of elements and material indices')
        N_node = coordinates.shape[0]
        elem_type = ELEM_TYPES.get((dim := coordinates.shape[1],
                                    nodes_per_elem := connectivity.shape[1]))
        if elem_type is None:
            raise ValueError(f'unsupported (dim, nodes_per_elem) = ({dim}, {nodes_per_elem})')

        self.id = str(uuid.uuid4())
        gmsh.model.add(self.id)
        gmsh.model.set_current(self.id)

        # gmsh always expects 3D coordinates, pad 2D with zeros
        coords_3d = np.pad(coordinates, ((0, 0), (0, 3 - coordinates.shape[1])))

        # classify each node on the entity of the first material using it, so that
        # no artificial node-holding entity with degenerate bounding box is needed;
        # nodes not referenced by any element end up on a separate holding entity
        owner = np.full(N_node, -1, dtype=np.int64)
        etags = []
        for i, mat_id in enumerate(np.unique(material)):
            etags.append(gmsh.model.add_discrete_entity(dim))
            used = np.unique(connectivity[material == mat_id])
            owner[used[owner[used] == -1]] = i

        if np.any(owner == -1):
            raise ValueError('dangling nodes')

        for i, mat_id in enumerate(np.unique(material)):
            sel = owner == i
            gmsh.model.mesh.add_nodes(dim, etags[i],
                                      np.flatnonzero(sel).astype(np.int64) + 1,
                                      coords_3d[sel].flatten())
            mask = material == mat_id
            gmsh.model.mesh.add_elements(dim, etags[i], [elem_type],
                                         [np.where(mask)[0].astype(np.int64) + 1],
                                         [connectivity[mask].flatten() + 1])
            gmsh.model.add_physical_group(dim, [etags[i]], tag=int(mat_id + 1))


    def __repr__(self) -> str:
        """
        Return repr(self).

        Give short, human-readable summary.
        """
        material = self.material
        mat_min = np.nanmin(material)
        mat_max = np.nanmax(material)
        mat_N   = np.unique(material).size
        return util.srepr([
               f'elements: {self.N_elements}',
               f'nodes:    {self.N_nodes}',
               f'# materials: {mat_N}' + ('' if mat_min == 0 and mat_max == mat_N-1 else
                                          f' (min: {mat_min}, max: {mat_max})')])


    def __del__(self):
        """Clean gmsh global state."""
        gmsh.model.set_current(self.id)
        gmsh.model.remove()


    @classmethod
    def from_grid(cls, grid) -> 'GeomMesh':
        """
        Create from a GeomGrid.

        Builds a hexahedral mesh and computes named boundary physical
        groups for all 6 faces, 12 edges, and 8 corners of the bounding box.

        Parameters
        ----------
        grid : damask.GeomGrid
            Source grid geometry.

        Returns
        -------
        mesh : damask.GeomMesh
            Mesh geometry with named boundary entities.
        """
        def _add_boundary(dim, name, tags):
            if dim == 0:  # corners become real points, otherwise $Entities shows fake origin coordinates
                etag = gmsh.model.geo.add_point(*coordinates[tags[0][0] - 1])
                gmsh.model.geo.synchronize()
            else:
                etag = gmsh.model.add_discrete_entity(dim)
            # element tags and physical group tags must be unique across the entire mesh
            elem_offset = max((int(b.max()) for b in gmsh.model.mesh.get_elements()[1] if len(b)), default=0)
            pg_offset = max((t for _, t in gmsh.model.get_physical_groups()), default=0)
            elem_tags = np.arange(elem_offset+1, elem_offset+len(tags)+1, dtype=np.int64)
            N_node = len(tags[0])
            etype = BOUNDARY_TYPES.get((dim, N_node), 1)
            gmsh.model.mesh.add_elements(dim, etag, [etype], [elem_tags],
                                         [np.array(tags, dtype=np.int64).flatten()])
            gmsh.model.add_physical_group(dim, [etag], tag=pg_offset+1, name=name)

        gmsh.option.set_number('General.Terminal', 0)

        coordinates = grid_.ravel(grid_.coordinates0_node(grid.cells, grid.size, grid.origin),
                                  flatten=True)
        node_idx = grid_.ravel_index(np.stack(np.meshgrid(*[np.arange(grid.cells[i] + 1) for i in range(3)],
                                                          indexing='ij'),
                                              axis=-1))
        # hex connectivity (0-based)
        connectivity = np.stack([node_idx[ :-1,  :-1,  :-1],
                                 node_idx[1:,    :-1,  :-1],
                                 node_idx[1:,   1:,    :-1],
                                 node_idx[ :-1, 1:,    :-1],
                                 node_idx[ :-1,  :-1, 1:],
                                 node_idx[1:,    :-1, 1:],
                                 node_idx[1:,   1:,   1:],
                                 node_idx[ :-1, 1:,   1:],
                                ], axis=-1).reshape(-1, 8, order='F')

        mesh = cls(coordinates, connectivity, grid.material.ravel(order='F'))

        # entities for boundary conditions
        gmsh.model.set_current(mesh.id)

        axis_names = ['x', 'y', 'z']
        for n_fixed in range(1, 4):
            dim = 3 - n_fixed
            for fixed_axes in itertools.combinations(range(3), n_fixed):
                for low in itertools.product([0, 1], repeat=n_fixed):
                    name = ''.join(f"{'-' if s == 0 else '+'}{axis_names[ax]}"
                                   for ax, s in zip(fixed_axes, low))
                    idx = [slice(None)] * 3
                    for ax, s in zip(fixed_axes, low):
                        idx[ax] = 0 if s == 0 else grid.cells[ax]                                   # type: ignore[assignment]

                    nodes = node_idx[tuple(idx)]
                    if dim == 2:  # face - quad elements
                        conn = np.stack([nodes[ :-1,  :-1],
                                         nodes[1:,    :-1],
                                         nodes[1:,   1:  ],
                                         nodes[ :-1, 1:  ]], axis=-1).reshape(-1, 4)

                    elif dim == 1:  # edge - line elements
                        conn = np.column_stack([nodes[:-1],
                                                nodes[1:]])

                    elif dim == 0:  # corner - point
                        conn = np.array([[nodes]])

                    _add_boundary(dim, name, conn+1)

        return mesh


    @classmethod
    def from_VTK(cls, vtk) -> 'GeomMesh':
        """
        Create from a damask.VTK.

        Parameters
        ----------
        vtk : damask.VTK
            VTK geometry of type UnstructuredGrid with cell data 'material'.

        Returns
        -------
        mesh : damask.GeomMesh
            Mesh geometry.
        """
        ug = vtk.vtk_data
        coordinates = vtk_to_numpy(ug.GetPoints().GetData())

        offsets = vtk_to_numpy(ug.GetCells().GetOffsetsArray())
        connectivity = vtk_to_numpy(ug.GetCells().GetConnectivityArray())
        N_cells = ug.GetNumberOfCells()
        cells = np.array([connectivity[offsets[i]:offsets[i + 1]]
                          for i in range(N_cells)], dtype=np.int64)

        # determine mesh dimension from cell type
        cell_type = ug.GetCellType(0)
        dim = {69: 2, 70: 2, 71: 3, 72: 3}.get(cell_type)                                          # Lagrange elements
        if dim is None:
            raise ValueError(f'unsupported VTK cell type {cell_type}')
        if cells.shape[1] != {69: 3, 70: 4, 71: 4, 72: 8}[cell_type]:
            raise ValueError('only linear elements are supported')

        return cls(coordinates[:, :dim], cells, vtk.get('material').astype(np.int64))


    def save(self, fname: Union[str, Path]) -> None:
        """
        Save mesh as gmsh file.

        Parameters
        ----------
        fname : str or pathlib.Path
            Output file name.
        """
        gmsh.model.set_current(self.id)
        gmsh.option.set_number('Mesh.SaveAll', 1)
        gmsh.write(str(Path(fname)))


    @staticmethod
    def load(fname: Union[str, Path]) -> 'GeomMesh':
        """
        Load mesh from gmsh file.

        Parameters
        ----------
        fname : str or pathlib.Path
            Input file name.

        Returns
        -------
        mesh : damask.GeomMesh
            Loaded mesh geometry.
        """
        mesh = object.__new__(GeomMesh)
        mesh.id = str(uuid.uuid4())
        gmsh.model.add(mesh.id)

        try:
            gmsh.merge(str(Path(fname)))

            dim = gmsh.model.get_dimension()
            if dim not in (2, 3):
                raise ValueError(f'unsupported mesh dimension {dim}')

            unique_types = np.unique(gmsh.model.mesh.get_elements(dim, -1)[0])
            if len(unique_types) != 1:
                raise ValueError(f'mesh contains more than one element type: {unique_types.tolist()}')
            if int(unique_types[0]) not in ELEM_TYPES.values():
                raise ValueError(f'unsupported element type {int(unique_types[0])}')

            return mesh
        finally:
            del mesh

    @property
    def coordinates(self) -> np.ndarray:
        """Node coordinates, shape (N, dim)."""
        gmsh.model.set_current(self.id)
        node_tags, coords, _ = gmsh.model.mesh.get_nodes()
        order = np.argsort(node_tags.astype(np.int64))

        return coords.reshape(-1, 3)[order]


    @property
    def connectivity(self) -> np.ndarray:
        """Element connectivity with 0-based node indices, shape (M, N_per_elem)."""
        gmsh.model.set_current(self.id)
        node_tags = gmsh.model.mesh.get_nodes()[0].astype(np.int64)
        order = np.argsort(node_tags)
        tag_to_idx = np.empty(node_tags.max() + 1, dtype=np.int64)
        tag_to_idx[node_tags[order]] = np.arange(len(node_tags))

        etypes, etags, enodes = gmsh.model.mesh.get_elements(gmsh.model.get_dimension(), -1)
        N_elem = len(etags[0])
        N_node = TYPE_TO_NODES.get(int(etypes[0]))

        return tag_to_idx[enodes[0].reshape(N_elem, N_node)]


    @property
    def material(self) -> np.ndarray:
        """Material ID per element (0-based), shape (M,)."""
        gmsh.model.set_current(self.id)
        dim = gmsh.model.get_dimension()
        result = []
        for d, etag in sorted(gmsh.model.get_entities(dim), key=lambda x: x[1]):
            pgs = gmsh.model.get_physical_groups_for_entity(d, int(etag))
            if len(pgs) == 0:
                continue
            mat_id = int(pgs[0]) - 1
            _, etags, _ = gmsh.model.mesh.get_elements(d, int(etag))
            n = len(etags[0])
            result.append(np.full(n, mat_id))

        return np.concatenate(result)

    @property
    def N_nodes(self) -> int:
        """Number of nodes."""
        gmsh.model.set_current(self.id)
        node_tags, _, _ = gmsh.model.mesh.get_nodes()
        return len(node_tags)

    @property
    def N_elements(self) -> int:
        """Number of elements."""
        gmsh.model.set_current(self.id)
        dim = gmsh.model.get_dimension()
        etypes, etags, _ = gmsh.model.mesh.get_elements(dim, -1)
        return sum(len(t) for t in etags)


    def show(self) -> None:
        """Visualize the mesh using gmsh."""
        gmsh.model.set_current(self.id)
        gmsh.fltk.run()
