import os

import pytest
import numpy as np

from damask import VTK

@pytest.fixture
def reference_dir(reference_dir_base):
    """Directory containing reference results."""
    return os.path.join(reference_dir_base,'Result')

class TestVTK:

    def test_rectilinearGrid(self,tmp_path):
        grid   = np.random.randint(5,10,3)*2
        size   = np.random.random(3) + 1.0
        origin = np.random.random(3)
        v = VTK.from_rectilinearGrid(grid,size,origin)
        s = v.__repr__()
        v.write(os.path.join(tmp_path,'rectilinearGrid'))
        v = VTK.from_file(os.path.join(tmp_path,'rectilinearGrid.vtr'))
        assert(s == v.__repr__())

    def test_polyData(self,tmp_path):
        points = np.random.rand(3,100)
        v = VTK.from_polyData(points)
        s = v.__repr__()
        v.write(os.path.join(tmp_path,'polyData'))
        v = VTK.from_file(os.path.join(tmp_path,'polyData.vtp'))
        assert(s == v.__repr__())

    @pytest.mark.parametrize('cell_type,n',[
                                            ('VTK_hexahedron',8),
                                            ('TETRA',4),
                                            ('quad',4),
                                            ('VTK_TRIANGLE',3)
                                            ]
                            )
    def test_unstructuredGrid(self,tmp_path,cell_type,n):
        nodes = np.random.rand(n,3)
        connectivity = np.random.choice(np.arange(n),n,False).reshape(-1,n)
        v = VTK.from_unstructuredGrid(nodes,connectivity,cell_type)
        s = v.__repr__()
        v.write(os.path.join(tmp_path,'unstructuredGrid'))
        v = VTK.from_file(os.path.join(tmp_path,'unstructuredGrid.vtu'))
        assert(s == v.__repr__())
