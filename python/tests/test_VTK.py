import os
import filecmp
import time
import string
import sys

import pytest
import numpy as np
import numpy.ma as ma
from vtkmodules.vtkCommonCore import vtkVersion

from damask import VTK
from damask import Table
from damask import Colormap

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'VTK'

@pytest.fixture
def default():
    """Simple VTK."""
    cells = np.array([5,6,7],int)
    size  = np.array([.6,1.,.5])
    return VTK.from_image_data(cells,size)

@pytest.fixture(autouse=True)
def _patch_execution_stamp(patch_execution_stamp):
    print('patched damask.util.execution_stamp')


@pytest.mark.parametrize('cmap',[Colormap.from_predefined('cividis'),'strain'])
@pytest.mark.skipif(sys.platform == 'win32', reason='DISPLAY has no effect on Windows OS')
def test_show(default,cmap,monkeypatch):
    monkeypatch.delenv('DISPLAY',raising=False)
    default.show(colormap=cmap)

def test_imageData(np_rng,tmp_path):
    cells = np_rng.integers(5,10,3)
    size = np_rng.random(3) + 0.1
    origin = np_rng.random(3) - 0.5
    v = VTK.from_image_data(cells,size,origin)
    string = str(v)
    string = v.as_ASCII()
    v.save(tmp_path/'imageData',False)
    vtr = VTK.load(tmp_path/'imageData.vti')
    with open(tmp_path/'imageData.vtk','w') as f:
        f.write(string)
    vtk = VTK.load(tmp_path/'imageData.vtk','VTK_imageData')
    assert (string == vtr.as_ASCII() == vtk.as_ASCII())

def test_rectilinearGrid(np_rng,tmp_path):
    grid = np.sort(np_rng.random((3,10)))
    v = VTK.from_rectilinear_grid(grid)
    string = str(v)
    string = v.as_ASCII()
    v.save(tmp_path/'rectilinearGrid',False)
    vtr = VTK.load(tmp_path/'rectilinearGrid.vtr')
    with open(tmp_path/'rectilinearGrid.vtk','w') as f:
        f.write(string)
    vtk = VTK.load(tmp_path/'rectilinearGrid.vtk','VTK_rectilinearGrid')
    assert (string == vtr.as_ASCII() == vtk.as_ASCII())

def test_polyData(np_rng,tmp_path):
    points = np_rng.random((100,3))
    v = VTK.from_poly_data(points)
    string = str(v)
    string = v.as_ASCII()
    v.save(tmp_path/'polyData',False)
    vtp = VTK.load(tmp_path/'polyData.vtp')
    with open(tmp_path/'polyData.vtk','w') as f:
        f.write(string)
    vtk = VTK.load(tmp_path/'polyData.vtk','polyData')
    assert (string == vtp.as_ASCII() == vtk.as_ASCII())

@pytest.mark.xfail(np.lib.NumpyVersion(vtkVersion.GetVTKVersion()) < '9.4.0',
                    reason = 'not available in VTK < 9.4')
def test_polyData_VTKHDF(np_rng,tmp_path):
    points = np_rng.random((100,3))
    v = VTK.from_poly_data(points)
    v.save_VTKHDF(tmp_path/'polyData')
    vtkhdf = VTK.load(tmp_path/'polyData.vtkhdf')
    assert (v.as_ASCII() == vtkhdf.as_ASCII())

# https://defelement.org/elements/lagrange.html
@pytest.mark.parametrize('cell_type,n_nodes',[
                                              ('VTK_hexahedron',lambda k: (k+1)**3),
                                              ('TETRA',lambda k: ((k+1)*(k+2)*(k+3))//6),
                                              ('quad',lambda k: (k+1)**2),
                                              ('VTK_TRIANGLE',lambda k: ((k+1)*(k+2))//2),
                                             ]
                        )
@pytest.mark.parametrize('order',[1,2,3,4,5])
def test_unstructuredGrid(np_rng,tmp_path,cell_type,n_nodes,order):
    n = n_nodes(order)
    nodes = np_rng.random((n,3))
    connectivity = np_rng.choice(np.arange(n),n,False).reshape(-1,n)
    v = VTK.from_unstructured_grid(nodes,connectivity,cell_type)
    string = str(v)
    string = v.as_ASCII()
    v.save(tmp_path/'unstructuredGrid',False)
    vtu = VTK.load(tmp_path/'unstructuredGrid.vtu')
    with open(tmp_path/'unstructuredGrid.vtk','w') as f:
        f.write(string)
    vtk = VTK.load(tmp_path/'unstructuredGrid.vtk','unstructuredgrid')
    assert (string == vtu.as_ASCII() == vtk.as_ASCII())

@pytest.mark.xfail(np.lib.NumpyVersion(vtkVersion.GetVTKVersion()) < '9.4.0',
                   reason = 'not available in VTK < 9.4')
@pytest.mark.parametrize('cell_type,n_nodes',[
                                              ('VTK_lagrange_hexa',lambda k: (k+1)**3),
                                              ('VTK_LAGRANGE_TETRAHEDRON',lambda k: ((k+1)*(k+2)*(k+3))//6),
                                              ('quadrilateral',lambda k: (k+1)**2),
                                              ('TRIANGLE',lambda k: ((k+1)*(k+2))//2),
                                             ]
                        )
@pytest.mark.parametrize('order',[1,2,3,4,5])
def test_unstructuredGrid_VTKHDF(np_rng,tmp_path,cell_type,n_nodes,order):
    n = n_nodes(order)
    nodes = np_rng.random((n,3))
    connectivity = np_rng.choice(np.arange(n),n,False).reshape(-1,n)
    v = VTK.from_unstructured_grid(nodes,connectivity,cell_type)
    v.save_VTKHDF(tmp_path/'unstructuredGrid')
    vtkhdf = VTK.load(tmp_path/'unstructuredGrid.vtkhdf')
    assert (v.as_ASCII() == vtkhdf.as_ASCII())


def test_parallel_out(np_rng,tmp_path):
    points = np_rng.random((102,3))
    v = VTK.from_poly_data(points)
    fname_s = tmp_path/'single.vtp'
    fname_p = tmp_path/'parallel.vtp'
    v.save(fname_s,False)
    v.save(fname_p,True)
    for i in range(10):
        if os.path.isfile(fname_p) and filecmp.cmp(fname_s,fname_p):
            assert True
            return
        time.sleep(.5)
    assert False

def test_compress(np_rng,tmp_path):
    points = np_rng.random((102,3))
    v = VTK.from_poly_data(points)
    fname_c = tmp_path/'compressed.vtp'
    fname_p = tmp_path/'plain.vtp'
    v.save(fname_c,parallel=False,compress=False)
    v.save(fname_p,parallel=False,compress=True)
    assert (VTK.load(fname_c).as_ASCII() == VTK.load(fname_p).as_ASCII())


@pytest.mark.parametrize('fname',['a','a.vtp','a.b','a.b.vtp'])
def test_filename_variations(np_rng,tmp_path,fname):
    points = np_rng.random((102,3))
    v = VTK.from_poly_data(points)
    v.save(tmp_path/fname)

@pytest.mark.parametrize('fname,dataset_type',[('a_file.vtk', None),
                                               ('a_file.vtk','vtk'),
                                               ('a_file.vtx', None)])
def test_invalid_dataset_type(tmp_path,fname,dataset_type):
    open(tmp_path/fname,'a').close()
    with pytest.raises(TypeError):
        VTK.load(tmp_path/fname,dataset_type)

def test_file_not_found():
    with pytest.raises(FileNotFoundError):
        VTK.load('/dev/null')

def test_add_extension(tmp_path,default):
    default.save(tmp_path/'default.txt',parallel=False)
    assert os.path.isfile(tmp_path/'default.txt.vti')


def test_invalid_get(default):
    with pytest.raises(KeyError):
        default.get('does_not_exist')

def test_invalid_set_shape(default):
    with pytest.raises(ValueError):
        default.set('valid',np.ones(3))

def test_invalid_set_missing_label(np_rng,default):
    data = np_rng.integers(9,size=np.prod(np.array(default.vtk_data.GetDimensions())-1))
    with pytest.raises(ValueError):
        default.set(data=data)

def test_invalid_set_type(default):
    with pytest.raises(TypeError):
        default.set(label='valid',data='invalid_type')
    with pytest.raises(TypeError):
        default.set(label='valid',table='invalid_type')

def test_invalid_set_dual(default):
    with pytest.raises(KeyError):
        default.set(label='valid',data=0,table=0)

@pytest.mark.parametrize('data_type,shape',[(float,(3,)),
                                            (float,(3,3)),
                                            (float,(1,)),
                                            (int,(4,)),
                                            (str,(1,))])
@pytest.mark.parametrize('N_values',[5*6*7,6*7*8])
def test_set_get(np_rng,default,data_type,shape,N_values):
    data = np.squeeze(np_rng.integers(0,100,(N_values,)+shape)).astype(data_type)
    new = default.set('data',data)
    assert (np.squeeze(data.reshape(N_values,-1)) == new.get('data')).all()


@pytest.mark.parametrize('shapes',[{'scalar':(1,),'vector':(3,),'tensor':(3,3)},
                                   {'vector':(6,),'tensor':(3,3)},
                                   {'tensor':(3,3),'scalar':(1,)}])
def test_set_table(np_rng,default,shapes):
    N = np_rng.choice([default.N_points,default.N_cells])
    d = dict()
    for k,s in shapes.items():
        d[k] = dict(shape = s,
                    data = np_rng.random(N*np.prod(s)).reshape((N,-1)))
    new = default.set(table=Table(shapes,np.column_stack([d[k]['data'] for k in shapes.keys()])))
    for k,s in shapes.items():
        assert np.allclose(np.squeeze(d[k]['data']),new.get(k),rtol=1e-7)


def test_set_masked(np_rng,default):
    data = np_rng.random((5*6*7,3))
    masked = ma.MaskedArray(data,mask=data<.4,fill_value=42.)
    mask_auto = default.set('D',masked)
    mask_manual = default.set('D',np.where(masked.mask,masked.fill_value,masked))
    assert mask_manual == mask_auto

@pytest.mark.parametrize('mode',['cells','points'])
def test_delete(np_rng,default,mode):
    data = np_rng.random(default.N_cells if mode == 'cells' else default.N_points).astype(np.float32)
    v = default.set('D',data)
    assert (v.get('D') == data).all()
    assert v.delete('D') == default

@pytest.mark.parametrize('data_type,shape',[(float,(3,)),
                                            (float,(3,3)),
                                            (float,(1,)),
                                            (int,(4,)),
                                            (str,(1,))])
@pytest.mark.parametrize('N_values',[5*6*7,6*7*8])
def test_labels(np_rng,default,data_type,shape,N_values):
    data = np.squeeze(np_rng.integers(0,100,(N_values,)+shape)).astype(data_type)
    ALPHABET = np.array(list(string.ascii_lowercase + ' '))
    label = ''.join(np_rng.choice(ALPHABET, size=10))
    new = default.set(label,data)
    if N_values == default.N_points: assert label in new.labels['Point Data']
    if N_values == default.N_cells:  assert label in new.labels['Cell Data']


def test_comments(tmp_path,default):
    default.comments += ['this is a comment']
    default.save(tmp_path/'with_comments',parallel=False)
    new = VTK.load(tmp_path/'with_comments.vti')
    assert new.comments == ['this is a comment']

def test_compare_reference_polyData(update,res_path):
    points=np.dstack((np.linspace(0.,1.,10),np.linspace(0.,2.,10),np.linspace(-1.,1.,10))).squeeze()
    polyData = VTK.from_poly_data(points).set('coordinates',points)
    if update:
        polyData.save(res_path/'polyData')
    else:
        reference = VTK.load(res_path/'polyData.vtp')
        assert polyData.as_ASCII() == reference.as_ASCII() and \
                np.allclose(polyData.get('coordinates'),points)

def test_compare_reference_rectilinearGrid(update,res_path):
    grid = [np.arange(4)**2.,
            np.arange(5)**2.,
            np.arange(6)**2.]                                               # ParaView renders tetrahedral meshing unless using float coordinates!
    coords = np.stack(np.meshgrid(*grid,indexing='ij'),axis=-1)
    c = coords[:-1,:-1,:-1,:].reshape(-1,3,order='F')
    n = coords[:,:,:,:].reshape(-1,3,order='F')
    rectilinearGrid = VTK.from_rectilinear_grid(grid) \
                    .set('cell',np.ascontiguousarray(c)) \
                    .set('node',np.ascontiguousarray(n))
    if update:
        rectilinearGrid.save(res_path/'rectilinearGrid')
    else:
        reference = VTK.load(res_path/'rectilinearGrid.vtr')
        assert rectilinearGrid.as_ASCII() == reference.as_ASCII() and \
                np.allclose(rectilinearGrid.get('cell'),c)
