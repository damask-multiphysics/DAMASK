import os
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

@pytest.fixture
def default():
    """Simple VTK."""
    grid = np.array([5,6,7],int)
    size = np.array([.6,1.,.5])
    return VTK.from_rectilinearGrid(grid,size)

class TestVTK:

    @pytest.fixture(autouse=True)
    def _execution_stamp(self, execution_stamp):
        print('patched damask.util.execution_stamp')

    def test_rectilinearGrid(self,tmp_path):
        grid   = np.random.randint(5,10,3)*2
        size   = np.random.random(3) + 1.0
        origin = np.random.random(3)
        v = VTK.from_rectilinearGrid(grid,size,origin)
        string = v.__repr__()
        v.write(tmp_path/'rectilinearGrid',False)
        vtr = VTK.from_file(tmp_path/'rectilinearGrid.vtr')
        with open(tmp_path/'rectilinearGrid.vtk','w') as f:
            f.write(string)
        vtk = VTK.from_file(tmp_path/'rectilinearGrid.vtk','VTK_rectilinearGrid')
        assert(string == vtr.__repr__() == vtk.__repr__())

    def test_polyData(self,tmp_path):
        points = np.random.rand(100,3)
        v = VTK.from_polyData(points)
        string = v.__repr__()
        v.write(tmp_path/'polyData',False)
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
        v.write(tmp_path/'unstructuredGrid',False)
        vtu = VTK.from_file(tmp_path/'unstructuredGrid.vtu')
        with open(tmp_path/'unstructuredGrid.vtk','w') as f:
            f.write(string)
        vtk = VTK.from_file(tmp_path/'unstructuredGrid.vtk','unstructuredgrid')
        assert(string == vtu.__repr__() == vtk.__repr__())


    def test_parallel_out(self,tmp_path):
        points = np.random.rand(102,3)
        v = VTK.from_polyData(points)
        fname_s = tmp_path/'single.vtp'
        fname_p = tmp_path/'parallel.vtp'
        v.write(fname_s,False)
        v.write(fname_p,True)
        for i in range(10):
            if os.path.isfile(fname_p) and filecmp.cmp(fname_s,fname_p):
                assert(True)
                return
            time.sleep(.5)
        assert(False)


    @pytest.mark.parametrize('name,dataset_type',[('this_file_does_not_exist.vtk', None),
                                                  ('this_file_does_not_exist.vtk','vtk'),
                                                  ('this_file_does_not_exist.vtx', None)])
    def test_invalid_dataset_type(self,name,dataset_type):
        with pytest.raises(TypeError):
            VTK.from_file(name,dataset_type)

    def test_invalid_extension_write(self,default):
        with pytest.raises(ValueError):
            default.write('default.txt')

    def test_invalid_get(self,default):
        with pytest.raises(ValueError):
            default.get('does_not_exist')

    @pytest.mark.parametrize('data,label',[(np.ones(3),'valid'),
                                           (np.ones(3),None)])
    def test_invalid_add(self,default,data,label):
        with pytest.raises(ValueError):
            default.add(np.ones(3),label)

    def test_invalid_add_type(self,default):
        with pytest.raises(TypeError):
            default.add('invalid_type','label')

    def test_comments(self,tmp_path,default):
        default.add_comments(['this is a comment'])
        default.write(tmp_path/'with_comments',parallel=False)
        new = VTK.from_file(tmp_path/'with_comments.vtr')
        assert new.get_comments() == ['this is a comment']


    def test_compare_reference_polyData(self,update,reference_dir,tmp_path):
        points=np.dstack((np.linspace(0.,1.,10),np.linspace(0.,2.,10),np.linspace(-1.,1.,10))).squeeze()
        polyData = VTK.from_polyData(points)
        polyData.add(points,'coordinates')
        if update:
             polyData.write(reference_dir/'polyData')
        else:
             reference = VTK.from_file(reference_dir/'polyData.vtp')
             assert polyData.__repr__() == reference.__repr__() and \
                    np.allclose(polyData.get('coordinates'),points)

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
             reference = VTK.from_file(reference_dir/'rectilinearGrid.vtr')
             assert rectilinearGrid.__repr__() == reference.__repr__() and \
                    np.allclose(rectilinearGrid.get('cell'),c)
