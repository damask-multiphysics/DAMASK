import os
import filecmp
import time

import pytest
import numpy as np
import numpy.ma as ma

from damask import VTK
from damask import grid_filters

@pytest.fixture
def ref_path(ref_path_base):
    """Directory containing reference results."""
    return ref_path_base/'VTK'

@pytest.fixture
def default():
    """Simple VTK."""
    cells = np.array([5,6,7],int)
    size  = np.array([.6,1.,.5])
    return VTK.from_rectilinear_grid(cells,size)

class TestVTK:

    @pytest.fixture(autouse=True)
    def _patch_execution_stamp(self, patch_execution_stamp):
        print('patched damask.util.execution_stamp')

    def test_rectilinearGrid(self,tmp_path):
        cells  = np.random.randint(5,10,3)*2
        size   = np.random.random(3) + 1.0
        origin = np.random.random(3)
        v = VTK.from_rectilinear_grid(cells,size,origin)
        string = v.__repr__()
        v.save(tmp_path/'rectilinearGrid',False)
        vtr = VTK.load(tmp_path/'rectilinearGrid.vtr')
        with open(tmp_path/'rectilinearGrid.vtk','w') as f:
            f.write(string)
        vtk = VTK.load(tmp_path/'rectilinearGrid.vtk','VTK_rectilinearGrid')
        assert(string == vtr.__repr__() == vtk.__repr__())

    def test_polyData(self,tmp_path):
        points = np.random.rand(100,3)
        v = VTK.from_poly_data(points)
        string = v.__repr__()
        v.save(tmp_path/'polyData',False)
        vtp = VTK.load(tmp_path/'polyData.vtp')
        with open(tmp_path/'polyData.vtk','w') as f:
            f.write(string)
        vtk = VTK.load(tmp_path/'polyData.vtk','polyData')
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
        v = VTK.from_unstructured_grid(nodes,connectivity,cell_type)
        string = v.__repr__()
        v.save(tmp_path/'unstructuredGrid',False)
        vtu = VTK.load(tmp_path/'unstructuredGrid.vtu')
        with open(tmp_path/'unstructuredGrid.vtk','w') as f:
            f.write(string)
        vtk = VTK.load(tmp_path/'unstructuredGrid.vtk','unstructuredgrid')
        assert(string == vtu.__repr__() == vtk.__repr__())


    def test_parallel_out(self,tmp_path):
        points = np.random.rand(102,3)
        v = VTK.from_poly_data(points)
        fname_s = tmp_path/'single.vtp'
        fname_p = tmp_path/'parallel.vtp'
        v.save(fname_s,False)
        v.save(fname_p,True)
        for i in range(10):
            if os.path.isfile(fname_p) and filecmp.cmp(fname_s,fname_p):
                assert(True)
                return
            time.sleep(.5)
        assert(False)

    def test_compress(self,tmp_path):
        points = np.random.rand(102,3)
        v = VTK.from_poly_data(points)
        fname_c = tmp_path/'compressed.vtp'
        fname_p = tmp_path/'plain.vtp'
        v.save(fname_c,parallel=False,compress=False)
        v.save(fname_p,parallel=False,compress=True)
        assert(VTK.load(fname_c).__repr__() == VTK.load(fname_p).__repr__())


    @pytest.mark.parametrize('fname',['a','a.vtp','a.b','a.b.vtp'])
    def test_filename_variations(self,tmp_path,fname):
        points = np.random.rand(102,3)
        v = VTK.from_poly_data(points)
        v.save(tmp_path/fname)

    @pytest.mark.parametrize('fname,dataset_type',[('a_file.vtk', None),
                                                   ('a_file.vtk','vtk'),
                                                   ('a_file.vtx', None)])
    def test_invalid_dataset_type(self,tmp_path,fname,dataset_type):
        open(tmp_path/fname,'a').close()
        with pytest.raises(TypeError):
            VTK.load(tmp_path/fname,dataset_type)

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            VTK.load('/dev/null')

    def test_add_extension(self,tmp_path,default):
        default.save(tmp_path/'default.txt',parallel=False)
        assert os.path.isfile(tmp_path/'default.txt.vtr')


    def test_invalid_get(self,default):
        with pytest.raises(ValueError):
            default.get('does_not_exist')

    def test_invalid_add_shape(self,default):
        with pytest.raises(ValueError):
            default.add(np.ones(3),'valid')

    def test_invalid_add_missing_label(self,default):
        data = np.random.randint(9,size=np.prod(np.array(default.vtk_data.GetDimensions())-1))
        with pytest.raises(ValueError):
            default.add(data)

    def test_invalid_add_type(self,default):
        with pytest.raises(TypeError):
            default.add('invalid_type','valid')

    @pytest.mark.parametrize('data_type,shape',[(float,(3,)),
                                                (float,(3,3)),
                                                (float,(1,)),
                                                (int,(4,)),
                                                (str,(1,))])
    @pytest.mark.parametrize('N_values',[5*6*7,6*7*8])
    def test_add_get(self,default,data_type,shape,N_values):
        data = np.squeeze(np.random.randint(0,100,(N_values,)+shape)).astype(data_type)
        default.add(data,'data')
        assert (np.squeeze(data.reshape(N_values,-1)) == default.get('data')).all()


    def test_add_masked(self,default):
        data = np.random.rand(5*6*7,3)
        masked = ma.MaskedArray(data,mask=data<.4,fill_value=42.)
        default.add(masked,'D')
        result_masked = str(default)
        default.add(np.where(masked.mask,masked.fill_value,masked),'D')
        assert result_masked == str(default)


    def test_comments(self,tmp_path,default):
        default.add_comments(['this is a comment'])
        default.save(tmp_path/'with_comments',parallel=False)
        new = VTK.load(tmp_path/'with_comments.vtr')
        assert new.get_comments() == ['this is a comment']

    def test_compare_reference_polyData(self,update,ref_path,tmp_path):
        points=np.dstack((np.linspace(0.,1.,10),np.linspace(0.,2.,10),np.linspace(-1.,1.,10))).squeeze()
        polyData = VTK.from_poly_data(points)
        polyData.add(points,'coordinates')
        if update:
             polyData.save(ref_path/'polyData')
        else:
             reference = VTK.load(ref_path/'polyData.vtp')
             assert polyData.__repr__() == reference.__repr__() and \
                    np.allclose(polyData.get('coordinates'),points)

    def test_compare_reference_rectilinearGrid(self,update,ref_path,tmp_path):
        cells = np.array([5,6,7],int)
        size  = np.array([.6,1.,.5])
        rectilinearGrid = VTK.from_rectilinear_grid(cells,size)
        c = grid_filters.coordinates0_point(cells,size).reshape(-1,3,order='F')
        n = grid_filters.coordinates0_node(cells,size).reshape(-1,3,order='F')
        rectilinearGrid.add(c,'cell')
        rectilinearGrid.add(n,'node')
        if update:
             rectilinearGrid.save(ref_path/'rectilinearGrid')
        else:
             reference = VTK.load(ref_path/'rectilinearGrid.vtr')
             assert rectilinearGrid.__repr__() == reference.__repr__() and \
                    np.allclose(rectilinearGrid.get('cell'),c)
