import numpy as np

from damask import Colormap

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
        rgbs = np.vstack((specials,np.random.rand(10000,3)))
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
