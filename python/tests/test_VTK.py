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
        string = v.__repr__()
        v.write(os.path.join(tmp_path,'rectilinearGrid'))
        vtr = VTK.from_file(os.path.join(tmp_path,'rectilinearGrid.vtr'))
        with open(os.path.join(tmp_path,'rectilinearGrid.vtk'),'w') as f:
            f.write(string)
        vtk = VTK.from_file(os.path.join(tmp_path,'rectilinearGrid.vtk'),'VTK_rectilinearGrid')
        assert(string == vtr.__repr__() == vtk.__repr__())

    def test_polyData(self,tmp_path):
        points = np.random.rand(3,100)
        v = VTK.from_polyData(points)
        string = v.__repr__()
        v.write(os.path.join(tmp_path,'polyData'))
        vtp = VTK.from_file(os.path.join(tmp_path,'polyData.vtp'))
        with open(os.path.join(tmp_path,'polyData.vtk'),'w') as f:
            f.write(string)
        vtk = VTK.from_file(os.path.join(tmp_path,'polyData.vtk'),'polyData')
        assert(string == vtp.__repr__() == vtk.__repr__())

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
        string = v.__repr__()
        v.write(os.path.join(tmp_path,'unstructuredGrid'))
        vtu = VTK.from_file(os.path.join(tmp_path,'unstructuredGrid.vtu'))
        with open(os.path.join(tmp_path,'unstructuredGrid.vtk'),'w') as f:
            f.write(string)
        vtk = VTK.from_file(os.path.join(tmp_path,'unstructuredGrid.vtk'),'unstructuredgrid')
        assert(string == vtu.__repr__() == vtk.__repr__())

    @pytest.mark.parametrize('name,dataset_type',[('this_file_does_not_exist.vtk',None),
                                                  ('this_file_does_not_exist.vtk','vtk'),
                                                  ('this_file_does_not_exist.vtx', None)])
    def test_invalid_dataset_type(self,dataset_type,name):
        with pytest.raises(TypeError):
            VTK.from_file('this_file_does_not_exist.vtk',dataset_type)
