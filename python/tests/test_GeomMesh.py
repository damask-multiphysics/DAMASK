# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np

from damask import GeomMesh
from damask import GeomGrid
from damask import VTK


# triangle (2D): 2 triangles sharing an edge
COORDINATES_TRIANGLE = np.array([
    [0.0, 0.0],
    [1.0, 0.0],
    [0.5, 1.0],
    [1.5, 1.0],
], dtype=np.float64)

CONNECTIVITY_TRIANGLE = np.array([
    [0, 1, 2],
    [1, 3, 2],
], dtype=np.int64)

MATERIAL_TRIANGLE = np.array([0, 1], dtype=np.int64)

CELL_TYPE_TRIANGLE = 'TRIANGLE'


# quadrilateral (2D): 2 quads sharing an edge
COORDINATES_QUAD = np.array([
    [0.0, 0.0],
    [1.0, 0.0],
    [2.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [2.0, 1.0],
], dtype=np.float64)

CONNECTIVITY_QUAD = np.array([
    [0, 1, 4, 3],
    [1, 2, 5, 4],
], dtype=np.int64)

MATERIAL_QUAD = np.array([0, 1], dtype=np.int64)

CELL_TYPE_QUAD = 'QUADRILATERAL'


# tetrahedron (3D): 2 tetrahedra from a cube split along one diagonal
COORDINATES_TET = np.array([
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [0.0, 1.0, 1.0],
    [1.0, 1.0, 1.0],
], dtype=np.float64)

CONNECTIVITY_TET = np.array([
    [0, 1, 3, 5],
    [0, 3, 2, 6],
    [0, 5, 4, 6],
    [5, 3, 7, 6],
    [0, 3, 5, 6],
], dtype=np.int64)

MATERIAL_TET = np.array([0, 0, 1, 1, 2], dtype=np.int64)

CELL_TYPE_TET = 'TETRAHEDRON'


# hexahedron (3D): 2 hexes sharing a face
COORDINATES_HEX = np.array([
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [2.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [1.0, 1.0, 0.0],
    [2.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [2.0, 0.0, 1.0],
    [0.0, 1.0, 1.0],
    [1.0, 1.0, 1.0],
    [2.0, 1.0, 1.0],
], dtype=np.float64)

CONNECTIVITY_HEX = np.array([
    [0, 1, 4, 3, 6, 7, 10, 9],
    [1, 2, 5, 4, 7, 8, 11, 10],
], dtype=np.int64)

MATERIAL_HEX = np.array([0, 1], dtype=np.int64)

CELL_TYPE_HEX = 'HEXAHEDRON'


meshes = [(COORDINATES_TRIANGLE,  CONNECTIVITY_TRIANGLE,  MATERIAL_TRIANGLE,  CELL_TYPE_TRIANGLE),
          (COORDINATES_QUAD,      CONNECTIVITY_QUAD,      MATERIAL_QUAD,      CELL_TYPE_QUAD),
          (COORDINATES_TET,       CONNECTIVITY_TET,       MATERIAL_TET,       CELL_TYPE_TET),
          (COORDINATES_HEX,       CONNECTIVITY_HEX,       MATERIAL_HEX,       CELL_TYPE_HEX)]

@pytest.fixture(autouse=True)
def _patch_execution_stamp(patch_execution_stamp):
    print('patched damask.util.execution_stamp')


@pytest.fixture(autouse=True)
def _patch_datetime_now(patch_datetime_now):
    print('patched datetime.datetime.now')


@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'GeomMesh'


if GeomMesh is None:
    pytest.skip('gmsh not installed', allow_module_level=True)
else:
    import gmsh


@pytest.mark.parametrize('mesh', meshes)
def test_init(mesh,assert_allclose):
    N_model0 = len(gmsh.model.list())

    coords, conn, mat, _ = mesh
    m = GeomMesh(coords, conn, mat)

    assert len(gmsh.model.list()) == N_model0 + 1
    assert m.N_nodes == coords.shape[0]
    assert m.N_elements == conn.shape[0]
    assert np.array_equal(m.connectivity, conn)
    assert np.array_equal(m.material, mat)
    assert_allclose(m.coordinates[:, :coords.shape[1]], coords)

@pytest.mark.parametrize('mesh', meshes)
def test_save_load(mesh,tmp_path,assert_allclose):
    N_model0 = len(gmsh.model.list())

    coords, conn, mat, _ = mesh
    m = GeomMesh(coords, conn, mat)
    m.save(tmp_path/'m.msh')
    m2 = GeomMesh.load(tmp_path/'m.msh')

    assert len(gmsh.model.list()) == N_model0 + 2
    assert m.N_nodes == m2.N_nodes
    assert m.N_elements == m2.N_elements
    assert np.array_equal(m.connectivity, m2.connectivity)
    assert np.array_equal(m.material, m2.material)
    assert_allclose(m.coordinates, m2.coordinates)


@pytest.mark.parametrize('mesh', meshes)
def test_from_vtk_roundtrip(mesh):
    N_model0 = len(gmsh.model.list())

    coords, conn, mat, cell_type = mesh
    m_init = GeomMesh(coords, conn, mat)

    coords_3d = np.pad(coords, ((0, 0), (0, 3 - coords.shape[1])))
    v = VTK.from_unstructured_grid(coords_3d, conn, cell_type)
    v = v.set('material', mat)
    m_vtk = GeomMesh.from_VTK(v)

    assert len(gmsh.model.list()) == N_model0 + 2
    assert m_init.N_nodes == m_vtk.N_nodes
    assert m_init.N_elements == m_vtk.N_elements
    assert np.allclose(m_init.coordinates, m_vtk.coordinates)
    assert np.array_equal(m_init.material, m_vtk.material)


def test_from_grid_roundtrip(np_rng):
    N_model0 = len(gmsh.model.list())

    size = (np_rng.random(3)+1.)*1e-5
    origin = np_rng.random(3)*1e-3
    cells = np_rng.integers(10,20,3)
    material = np_rng.integers(0,40,size=cells)
    g = GeomGrid(material,size,origin)
    m = GeomMesh.from_grid(g)

    assert len(gmsh.model.list()) == N_model0 + 1
    assert np.array_equal(np.sort(material.flatten()), np.sort(m.material))
    assert cells.prod() == m.N_elements
    assert (cells+1).prod() == m.N_nodes

def test_from_grid_regular(np_rng,assert_allclose):
    N_model0 = len(gmsh.model.list())

    size = (np_rng.random(3)+1.)*1e-5
    origin = np_rng.random(3)*1e-3
    cells = np_rng.integers(10,20,3)
    material = np.arange(cells.prod()).reshape(cells)
    g = GeomGrid(material,size,origin)
    m = GeomMesh.from_grid(g)

    assert len(gmsh.model.list()) == N_model0 + 1
    assert np.array_equal(material.flatten(), m.material)
    assert cells.prod() == m.N_elements
    assert (cells+1).prod() == m.N_nodes

    N_material = cells.prod()
    groups = {}
    for d,t in gmsh.model.get_physical_groups():
        groups.setdefault(d,[]).append(t)
    assert sorted(groups[0]) == list(range(N_material+19,N_material+27))                            # 8 corners
    assert sorted(groups[1]) == list(range(N_material+7,N_material+19))                             # 12 edges
    assert sorted(groups[2]) == list(range(N_material+1,N_material+7))                              # 6 faces

    for d,t in gmsh.model.get_entities(0):                                                          # corners have true coordinates
        pg = int(gmsh.model.get_physical_groups_for_entity(d,t)[0])
        name = gmsh.model.get_physical_name(d,pg)
        corner = [origin[i] + (0 if name[2*i] == '-' else size[i]) for i in range(3)]
        assert_allclose(gmsh.model.get_bounding_box(d,t),corner+corner)

def test_element_tag_count(np_rng):
    cells = np_rng.integers(2,6,3)
    N_material = np_rng.integers(1,cells.prod()+1)
    material = np_rng.integers(0,N_material,cells)
    g = GeomGrid(material,np_rng.random(3)+0.1)
    m = GeomMesh.from_grid(g)
    gmsh.model.set_current(m.id)

    N_node = {1: 2, 3: 4, 5: 8, 15: 1}                                                              # per element type
    tags, elements = [], []
    for d,_ in gmsh.model.get_entities():
        types, elem_tags, node_tags = gmsh.model.mesh.get_elements(d,_)
        assert len(elem_tags) == len(types) == len(node_tags) == 1                                  # single block
        n = N_node[int(types[0])]
        assert node_tags[0].size == n*len(elem_tags[0])
        tags.extend(elem_tags[0].tolist())
        elements.extend((d,int(types[0]),tuple(sorted(e.tolist())))
                         for e in np.split(node_tags[0],len(elem_tags[0])))
    assert len(tags) == len(set(tags))                                                              # tags unique
    assert len(set(elements)) == len(elements)                                                      # no element twice
    assert len(elements) == cells.prod() \
                           + 2*(cells[0]*cells[1]+cells[0]*cells[2]+cells[1]*cells[2]) \
                           + 4*(cells.sum()) + 8                                                    # hex+quad+line+point
    assert sorted(tags) == list(range(1,len(elements)+1))                                           # consecutive

def test_from_vtk_quadratic():
    N_model0 = len(gmsh.model.list())

    coords = np.array([[x, y, z]
                       for x in (0., 0.5, 1.)
                       for y in (0., 0.5, 1.)
                       for z in (0., 0.5, 1.)], dtype=np.float64)
    connectivity = np.arange(27, dtype=np.int64).reshape(1, -1)
    v = VTK.from_unstructured_grid(coords, connectivity, 'HEXAHEDRON')
    v = v.set('material', np.array([0]))

    with pytest.raises(ValueError, match='linear'):
        GeomMesh.from_VTK(v)
    assert len(gmsh.model.list()) == N_model0

@pytest.mark.parametrize('fname', ['mixed_2D.msh', 'mixed_3D.msh'])
def test_load_mixed_element_types(res_path, fname):
    N_model0 = len(gmsh.model.list())

    with pytest.raises(ValueError, match='more than one element type'):
        GeomMesh.load(res_path/fname)
    assert len(gmsh.model.list()) == N_model0

def test_load_unsupported_element_type(res_path):
    N_model0 = len(gmsh.model.list())

    with pytest.raises(ValueError, match='unsupported element type'):
        GeomMesh.load(res_path/'prism.msh')
    assert len(gmsh.model.list()) == N_model0


def test_cleanup():
    N_model0 = len(gmsh.model.list())
    coords, conn, mat, _ = meshes[0]
    m1 = GeomMesh(coords, conn, mat)
    assert len(gmsh.model.list()) == N_model0+1
    m2 = GeomMesh(coords, conn, mat)
    assert len(gmsh.model.list()) == N_model0+2
    m1 = m2
    assert len(gmsh.model.list()) == N_model0+1
    del m1
    del m2
    assert len(gmsh.model.list()) == N_model0
