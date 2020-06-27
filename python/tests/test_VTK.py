import filecmp
import time

import pytest
import numpy as np

from damask import VTK
from damask import grid_filters

@pytest.fixture
def reference_dir(reference_dir_base):
    """Directory containing reference results."""
    return reference_dir_base/'VTK'

class TestVTK:

    def test_rectilinearGrid(self,tmp_path):
        grid   = np.random.randint(5,10,3)*2
        size   = np.random.random(3) + 1.0
        origin = np.random.random(3)
        v = VTK.from_rectilinearGrid(grid,size,origin)
        string = v.__repr__()
        v.write(tmp_path/'rectilinearGrid')
        vtr = VTK.from_file(tmp_path/'rectilinearGrid.vtr')
        with open(tmp_path/'rectilinearGrid.vtk','w') as f:
            f.write(string)
        vtk = VTK.from_file(tmp_path/'rectilinearGrid.vtk','VTK_rectilinearGrid')
        assert(string == vtr.__repr__() == vtk.__repr__())

    def test_polyData(self,tmp_path):
        points = np.random.rand(3,100)
        v = VTK.from_polyData(points)
        string = v.__repr__()
        v.write(tmp_path/'polyData')
        vtp = VTK.from_file(tmp_path/'polyData.vtp')
        with open(tmp_path/'polyData.vtk','w') as f:
            f.write(string)
        vtk = VTK.from_file(tmp_path/'polyData.vtk','polyData')
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
        v.write(tmp_path/'unstructuredGrid')
        vtu = VTK.from_file(tmp_path/'unstructuredGrid.vtu')
        with open(tmp_path/'unstructuredGrid.vtk','w') as f:
            f.write(string)
        vtk = VTK.from_file(tmp_path/'unstructuredGrid.vtk','unstructuredgrid')
        assert(string == vtu.__repr__() == vtk.__repr__())

    @pytest.mark.parametrize('name,dataset_type',[('this_file_does_not_exist.vtk',None),
                                                  ('this_file_does_not_exist.vtk','vtk'),
                                                  ('this_file_does_not_exist.vtx', None)])
    def test_invalid_dataset_type(self,dataset_type,name):
        with pytest.raises(TypeError):
            VTK.from_file('this_file_does_not_exist.vtk',dataset_type)


    def test_compare_reference_polyData(self,update,reference_dir,tmp_path):
        points=np.dstack((np.linspace(0.,1.,10),np.linspace(0.,2.,10),np.linspace(-1.,1.,10))).squeeze()
        polyData = VTK.from_polyData(points)
        polyData.add(points,'coordinates')
        if update:
             polyData.write(reference_dir/'polyData')
        else:
             polyData.write(tmp_path/'polyData')
             time.sleep(.5)
             assert filecmp.cmp(tmp_path/'polyData.vtp',reference_dir/'polyData.vtp')

    def test_compare_reference_rectilinearGrid(self,update,reference_dir,tmp_path):
        grid = np.array([5,6,7],int)
        size = np.array([.6,1.,.5])
        rectilinearGrid = VTK.from_rectilinearGrid(grid,size)
        c = grid_filters.cell_coord0(grid,size).reshape(-1,3,order='F')
        n = grid_filters.node_coord0(grid,size).reshape(-1,3,order='F')
        rectilinearGrid.add(c,'cell')
        rectilinearGrid.add(n,'node')
        if update:
             rectilinearGrid.write(reference_dir/'rectilinearGrid')
        else:
             rectilinearGrid.write(tmp_path/'rectilinearGrid')
             time.sleep(.5)
             assert filecmp.cmp(tmp_path/'rectilinearGrid.vtr',reference_dir/'rectilinearGrid.vtr')

