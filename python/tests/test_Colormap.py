import os
import filecmp
import time

import numpy as np
import pytest
from PIL import Image
from PIL import ImageChops

import damask
from damask import Colormap

@pytest.fixture
def reference_dir(reference_dir_base):
    """Directory containing reference results."""
    return reference_dir_base/'Colormap'

class TestColormap:

    def test_conversion(self):

        specials = np.array([[0.,0.,0.],
                             [1.,0.,0.],
                             [0.,1.,0.],
                             [0.,0.,1.],
                             [1.,1.,0.],
                             [0.,1.,1.],
                             [1.,0.,1.],
                             [1.,1.,1.]
                             ])
        rgbs = np.vstack((specials,np.random.rand(100,3)))
        pass # class not integrated
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


    @pytest.mark.parametrize('format',['ASCII','paraview','GOM','Gmsh'])
    @pytest.mark.parametrize('model',['rgb','hsv','hsl','xyz','lab','msh'])
    def test_from_range(self,model,format,tmpdir):
        N = np.random.randint(2,256)
        c = Colormap.from_range(np.random.rand(3),np.random.rand(3),model=model,N=N)
        c.to_file(tmpdir/'color_out',format=format)

    @pytest.mark.parametrize('format',['ASCII','paraview','GOM','Gmsh'])
    @pytest.mark.parametrize('name',['strain','gnuplot','Greys','PRGn','viridis'])
    def test_from_predefined(self,name,format,tmpdir):
        N = np.random.randint(2,256)
        c = Colormap.from_predefined(name,N)
        os.chdir(tmpdir)
        c.to_file(format=format)

    @pytest.mark.parametrize('format,name',[('ASCII','test.txt'),
                                            ('paraview','test.json'),
                                            ('GOM','test.legend'),
                                            ('Gmsh','test.msh')
                                           ])
    def test_write_filehandle(self,format,name,tmpdir):
        c = Colormap.from_predefined('Dark2')
        fname = tmpdir/name
        with open(fname,'w') as f:
            c.to_file(f,format=format)
        for i in range(10):
            if fname.exists(): return
            time.sleep(.5)
        assert False

    def test_write_invalid_format(self):
        c = Colormap.from_predefined('Dark2')
        with pytest.raises(ValueError):
            c.to_file(format='invalid')

    @pytest.mark.parametrize('model',['rgb','hsv','hsl','lab','invalid'])
    def test_invalid_color(self,model):
        with pytest.raises(ValueError):
            c = Colormap.from_range(-2.+np.random.rand(3),np.random.rand(3),N=10,model=model)      # noqa

    def test_reversed(self):
        c_1 = Colormap.from_predefined('stress')
        c_2 = c_1.reversed()
        assert (not np.allclose(c_1.colors,c_2.colors)) and \
                np.allclose(c_1.colors,c_2.reversed().colors)

    def test_invert(self):
        c_1 = Colormap.from_predefined('strain')
        c_2 = ~c_1
        assert (not np.allclose(c_1.colors,c_2.colors)) and \
                np.allclose(c_1.colors,(~c_2).colors)

    def test_add(self):
        c = Colormap.from_predefined('jet')
        c += c
        assert (np.allclose(c.colors[:len(c.colors)//2],c.colors[len(c.colors)//2:]))

    @pytest.mark.parametrize('bounds',[None,[2,10]])
    def test_shade(self,reference_dir,update,bounds):
        data = np.add(*np.indices((10, 11)))
        img_current = Colormap.from_predefined('orientation').shade(data,bounds=bounds)
        if update:
            img_current.save(reference_dir/f'shade_{bounds}.png')
        else:
            img_reference = Image.open(reference_dir/f'shade_{bounds}.png')
            diff = ImageChops.difference(img_reference.convert('RGB'),img_current.convert('RGB'))
            assert not diff.getbbox()

    def test_list(self):
        Colormap.list_predefined()

    @pytest.mark.parametrize('format,ext',[('ASCII','.txt'),
                                           ('paraview','.json'),
                                           ('GOM','.legend'),
                                           ('Gmsh','.msh')
                                          ])
    def test_compare_reference(self,format,ext,tmpdir,reference_dir,update,monkeypatch):
        monkeypatch.setattr(damask, 'version', pytest.dummy_version)
        name = 'binary'
        c = Colormap.from_predefined(name)
        if update:
            os.chdir(reference_dir)
            c.to_file(format=format)
        else:
            os.chdir(tmpdir)
            c.to_file(format=format)
            time.sleep(.5)
            assert filecmp.cmp(tmpdir/(name+ext),reference_dir/(name+ext))
