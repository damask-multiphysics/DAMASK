import os
import filecmp
import time

import numpy as np
import pytest
from PIL import Image
from PIL import ImageChops

from damask import Colormap

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'Colormap'

@pytest.fixture(autouse=True)
def _patch_execution_stamp(patch_execution_stamp):
    print('patched damask.util.execution_stamp')


def test_repr(patch_plt_show):
    print(Colormap.from_predefined('stress'))

def test_conversion(np_rng):
    specials = np.array([[0.,0.,0.],
                         [1.,0.,0.],
                         [0.,1.,0.],
                         [0.,0.,1.],
                         [1.,1.,0.],
                         [0.,1.,1.],
                         [1.,0.,1.],
                         [1.,1.,1.]
                        ])
    rgbs = np.vstack((specials,np_rng.random((100,3))))
    for rgb in rgbs:
        print('rgb',rgb)

        # rgb2hsv2rgb
        hsv = Colormap._rgb2hsv(rgb)
        print('hsv',hsv)
        assert np.allclose(Colormap._hsv2rgb(hsv),rgb)

        # rgb2hsl2rgb
        hsl = Colormap._rgb2hsl(rgb)
        print('hsl',hsl)
        assert np.allclose(Colormap._hsl2rgb(hsl),rgb)

        # rgb2xyz2rgb
        xyz = Colormap._rgb2xyz(rgb)
        print('xyz',xyz)
        assert np.allclose(Colormap._xyz2rgb(xyz),rgb,atol=1.e-6,rtol=0)

        # xyz2lab2xyz
        lab = Colormap._xyz2lab(xyz)
        print('lab',lab)
        assert np.allclose(Colormap._lab2xyz(lab),xyz)

        # lab2msh2lab
        msh = Colormap._lab2msh(lab)
        print('msh',msh)
        assert np.allclose(Colormap._msh2lab(msh),lab)

        # lab2rgb2lab
        assert np.allclose(Colormap._rgb2lab(Colormap._lab2rgb(lab)),lab,atol=1.e-6,rtol=0)

        # rgb2msh2rgb
        assert np.allclose(Colormap._msh2rgb(Colormap._rgb2msh(rgb)),rgb,atol=1.e-6,rtol=0)

        # hsv2msh
        assert np.allclose(Colormap._hsv2msh(hsv),msh,atol=1.e-6,rtol=0)

        # hsl2msh
        assert np.allclose(Colormap._hsv2msh(hsv),msh,atol=1.e-6,rtol=0)

        # xyz2msh
        assert np.allclose(Colormap._xyz2msh(xyz),msh,atol=1.e-6,rtol=0)

def test_eq():
    assert Colormap.from_predefined('strain') == Colormap.from_predefined('strain')
    assert Colormap.from_predefined('strain') != Colormap.from_predefined('stress')
    assert Colormap.from_predefined('strain',N=128) != Colormap.from_predefined('strain',N=64)
    assert not Colormap.from_predefined('strain',N=128) == 1

@pytest.mark.parametrize('low,high',[((0,0,0),(1,1,1)),
                                     ([0,0,0],[1,1,1])])
def test_from_range_types(low,high):
    assert Colormap.from_range(low,high) == Colormap.from_range(np.array(low),np.array(high))

@pytest.mark.parametrize('format',['ASCII','paraview','GOM','gmsh'])
@pytest.mark.parametrize('model',['rgb','hsv','hsl','xyz','lab','msh'])
def test_from_range(np_rng,model,format,tmp_path):
    N = np_rng.integers(2,256)
    c = Colormap.from_range(np_rng.random(3),np_rng.random(3),model=model,N=N)                      # noqa
    eval(f'c.save_{format}(tmp_path/"color_out")')

@pytest.mark.parametrize('format',['ASCII','paraview','GOM','gmsh'])
@pytest.mark.parametrize('name',['strain','gnuplot','Greys','PRGn','viridis'])
def test_from_predefined(np_rng,name,format,tmp_path):
    N = np_rng.integers(2,256)
    c = Colormap.from_predefined(name,N)                                                            # noqa
    os.chdir(tmp_path)
    eval(f'c.save_{format}()')

@pytest.mark.parametrize('format,name',[('ASCII','test.txt'),
                                        ('paraview','test.json'),
                                        ('GOM','test.legend'),
                                        ('gmsh','test.msh')
                                        ])
def test_write_filehandle(format,name,tmp_path):
    c = Colormap.from_predefined('Dark2')                                                           # noqa
    fname = tmp_path/name
    with open(fname,'w') as f:                                                                      # noqa
        eval(f'c.save_{format}(f)')
    for i in range(10):
        if fname.exists(): return
        time.sleep(.5)
    assert False

@pytest.mark.parametrize('model',['rgb','hsv','hsl','lab','invalid'])
def test_invalid_color(np_rng,model):
    with pytest.raises(ValueError):
        c = Colormap.from_range(-2.+np_rng.random(3),np_rng.random(3),N=10,model=model)             # noqa

def test_reversed():
    c_1 = Colormap.from_predefined('stress')
    c_2 = c_1.reversed()
    assert (not np.allclose(c_1.colors,c_2.colors)) and \
                np.allclose(c_1.colors,c_2.reversed().colors)

def test_invert():
    c_1 = Colormap.from_predefined('strain')
    c_2 = ~c_1
    assert (not np.allclose(c_1.colors,  c_2.colors)) and \
                np.allclose(c_1.colors,(~c_2).colors)

def test_add():
    c = Colormap.from_predefined('jet')
    c += c
    assert (np.allclose(c.colors[:len(c.colors)//2],c.colors[len(c.colors)//2:]))

def test_mul():
    c = o = Colormap.from_predefined('jet')
    o *= 2
    assert c+c == o

@pytest.mark.parametrize('N,cmap,at,result',[
        (8,'gray',0.5,[0.5,0.5,0.5]),
        (17,'gray',0.5,[0.5,0.5,0.5]),
        (17,'gray',[0.5,0.75],[[0.5,0.5,0.5],[0.75,0.75,0.75]]),
        ])
def test_at_value(N, cmap, at, result):
    assert np.allclose(Colormap.from_predefined(cmap,N=N).at(at)[...,:3],
                       result,
                       rtol=0.005)

@pytest.mark.parametrize('bounds',[None,[2,10]])
def test_shade(res_path,update,bounds):
    data = np.add(*np.indices((10, 11)))
    img_current = Colormap.from_predefined('orientation').shade(data,bounds=bounds)
    if update:
        img_current.save(res_path/f'shade_{bounds}.png')
    else:
        img_reference = Image.open(res_path/f'shade_{bounds}.png')
        diff = ImageChops.difference(img_reference.convert('RGB'),img_current.convert('RGB'))
        assert not diff.getbbox()

def test_predefined():
    assert (isinstance(Colormap.predefined,dict))

@pytest.mark.parametrize('format,ext',[('ASCII','.txt'),
                                       ('paraview','.json'),
                                       ('GOM','.legend'),
                                       ('gmsh','.msh')
                                       ])
def test_compare_reference(format,ext,tmp_path,res_path,update):
    name = 'binary'
    c = Colormap.from_predefined(name)                                                              # noqa
    if update:
        os.chdir(res_path)
        eval(f'c.save_{format}()')
    else:
        os.chdir(tmp_path)
        eval(f'c.save_{format}()')
        time.sleep(.5)
        assert filecmp.cmp(tmp_path/(name+ext),res_path/(name+ext))
