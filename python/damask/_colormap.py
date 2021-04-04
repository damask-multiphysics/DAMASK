import os
import json
import functools

import numpy as np
import matplotlib as mpl
if os.name == 'posix' and 'DISPLAY' not in os.environ:
    mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from PIL import Image

from . import util
from . import Table

_eps   = 216./24389.
_kappa = 24389./27.
_ref_white = np.array([.95047, 1.00000, 1.08883])                                                   # Observer = 2, Illuminant = D65

# ToDo (if needed)
# - support alpha channel (paraview/ASCII/input)
# - support NaN color (paraview)

class Colormap(mpl.colors.ListedColormap):
    """
    Enhance matplotlib colormap functionality to be used within DAMASK.

    References
    ----------
    K. Moreland, Proceedings of the 5th International Symposium on Advances in Visual Computing, 2009
    https://doi.org/10.1007/978-3-642-10520-3_9

    P. Eisenlohr et al., International Journal of Plasticity 46:37–53, 2013
    https://doi.org/10.1016/j.ijplas.2012.09.012

    Matplotlib colormaps overview
    https://matplotlib.org/tutorials/colors/colormaps.html

    """

    def __add__(self,other):
        """Concatenate."""
        return Colormap(np.vstack((self.colors,other.colors)),
                        f'{self.name}+{other.name}')

    def __iadd__(self,other):
        """Concatenate (in-place)."""
        return self.__add__(other)

    def __invert__(self):
        """Reverse."""
        return self.reversed()

    def __repr__(self):
        """Show as matplotlib figure."""
        fig = plt.figure(self.name,figsize=(5,.5))
        ax1 = fig.add_axes([0, 0, 1, 1])
        ax1.set_axis_off()
        ax1.imshow(np.linspace(0,1,self.N).reshape(1,-1),
                   aspect='auto', cmap=self, interpolation='nearest')
        plt.show(block = False)
        return 'Colormap: '+self.name


    @staticmethod
    def from_range(low,high,name='DAMASK colormap',N=256,model='rgb'):
        """
        Create a perceptually uniform colormap between given (inclusive) bounds.

        Colors are internally stored as R(ed) G(green) B(lue) values.
        The colormap can be used in matplotlib/seaborn or exported to
        file for external use.

        Parameters
        ----------
        low : numpy.ndarray of shape (3)
            Color definition for minimum value.
        high : numpy.ndarray of shape (3)
            Color definition for maximum value.
        N : int, optional
            The number of color quantization levels. Defaults to 256.
        name : str, optional
            The name of the colormap. Defaults to `DAMASK colormap`.
        model : {'rgb', 'hsv', 'hsl', 'xyz', 'lab', 'msh'}
            Colormodel used for input color definitions. Defaults to `rgb`.
            The available color models are:
            - 'rgb': R(ed) G(green) B(lue).
            - 'hsv': H(ue) S(aturation) V(alue).
            - 'hsl': H(ue) S(aturation) L(uminance).
            - 'xyz': CIE Xyz.
            - 'lab': CIE Lab.
            - 'msh': Msh (for perceptual uniform interpolation).

        """
        low_high = np.vstack((low,high))
        if   model.lower() == 'rgb':
            if np.any(low_high<0) or np.any(low_high>1):
                raise ValueError(f'RGB color {low} | {high} are out of range.')

            low_,high_ = map(Colormap._rgb2msh,low_high)

        elif model.lower() == 'hsv':
            if np.any(low_high<0) or np.any(low_high[:,1:3]>1) or np.any(low_high[:,0]>360):
                raise ValueError(f'HSV color {low} | {high} are out of range.')

            low_,high_ = map(Colormap._hsv2msh,low_high)

        elif model.lower() == 'hsl':
            if np.any(low_high<0) or np.any(low_high[:,1:3]>1) or np.any(low_high[:,0]>360):
                raise ValueError(f'HSL color {low} | {high} are out of range.')

            low_,high_ = map(Colormap._hsl2msh,low_high)

        elif model.lower() == 'xyz':

            low_,high_ = map(Colormap._xyz2msh,low_high)

        elif model.lower() == 'lab':
            if np.any(low_high[:,0]<0):
                raise ValueError(f'CIE Lab color {low} | {high} are out of range.')

            low_,high_ = map(Colormap._lab2msh,low_high)

        elif model.lower() == 'msh':
            low_,high_ = low_high[0],low_high[1]

        else:
            raise ValueError(f'Invalid color model: {model}.')

        msh = map(functools.partial(Colormap._interpolate_msh,low=low_,high=high_),np.linspace(0,1,N))
        rgb = np.array(list(map(Colormap._msh2rgb,msh)))

        return Colormap(rgb,name=name)


    @staticmethod
    def from_predefined(name,N=256):
        """
        Select from a set of predefined colormaps.

        Predefined colormaps include native matplotlib colormaps
        and common DAMASK colormaps.

        Parameters
        ----------
        name : str
            The name of the colormap.
        N : int, optional
           The number of color quantization levels. Defaults to 256.
           This parameter is not used for matplotlib colormaps
           that are of type `ListedColormap`.

        """
        # matplotlib presets
        try:
            colormap = cm.__dict__[name]
            return Colormap(np.array(list(map(colormap,np.linspace(0,1,N)))
                                     if isinstance(colormap,mpl.colors.LinearSegmentedColormap) else
                                     colormap.colors),
                            name=name)
        except KeyError:
            # DAMASK presets
            definition = Colormap._predefined_DAMASK[name]
            return Colormap.from_range(definition['low'],definition['high'],name,N)


    def shade(self,field,bounds=None,gap=None):
        """
        Generate PIL image of 2D field using colormap.

        Parameters
        ----------
        field : numpy.array of shape (:,:)
            Data to be shaded.
        bounds : iterable of len (2), optional
            Colormap value range (low,high).
        gap : field.dtype, optional
            Transparent value. NaN will always be rendered transparent.

        Returns
        -------
        PIL.Image
            RGBA image of shaded data.

        """
        N = len(self.colors)
        mask = np.logical_not(np.isnan(field) if gap is None else \
               np.logical_or (np.isnan(field), field == gap))                                       # mask NaN (and gap if present)

        lo,hi = (field[mask].min(),field[mask].max()) if bounds is None else \
                (min(bounds[:2]),max(bounds[:2]))

        delta,avg = hi-lo,0.5*(hi+lo)

        if delta * 1e8 <= avg:                                                                      # delta is similar to numerical noise
            hi,lo = hi+0.5*avg,lo-0.5*avg                                                           # extend range to have actual data centered within

        return Image.fromarray(
            (np.dstack((
                        self.colors[(np.round(np.clip((field-lo)/(hi-lo),0.0,1.0)*(N-1))).astype(np.uint16),:3],
                        mask.astype(float)
                       )
                      )*255
            ).astype(np.uint8),
            mode='RGBA')


    def reversed(self,name=None):
        """
        Reverse.

        Parameters
        ----------
        name : str, optional
            The name for the reversed colormap.
            A name of None will be replaced by the name of the parent colormap + "_r".

        Returns
        -------
        damask.Colormap
            The reversed colormap.

        """
        rev = super(Colormap,self).reversed(name)
        return Colormap(np.array(rev.colors),rev.name[:-4] if rev.name.endswith('_r_r') else rev.name)


    def _get_file_handle(self,fname,extension):
        """
        Provide file handle.

        Parameters
        ----------
        fname : file, str, pathlib.Path, or None
            Filename or filehandle, will be name of the colormap+extension if None.

        extension: str
            Extension of the filename.

        Returns
        -------
        f
            File handle

        """
        if fname is None:
            fhandle = open(self.name.replace(' ','_')+'.'+extension,'w',newline='\n')
        else:
            try:
                fhandle = open(fname,'w',newline='\n')
            except TypeError:
                fhandle = fname

        return fhandle


    def save_paraview(self,fname=None):
        """
        Save as JSON file for use in Paraview.

        Parameters
        ----------
        fname : file, str, or pathlib.Path, optional
            Filename to store results. If not given, the filename will
            consist of the name of the colormap and extension '.json'.

        """
        colors = []
        for i,c in enumerate(np.round(self.colors,6).tolist()):
            colors+=[i]+c

        out = [{
                'Creator':util.execution_stamp('Colormap'),
                'ColorSpace':'RGB',
                'Name':self.name,
                'DefaultMap':True,
                'RGBPoints':colors
               }]

        json.dump(out,self._get_file_handle(fname,'json'),indent=4)


    def save_ASCII(self,fname=None):
        """
        Save as ASCII file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path, optional
            Filename to store results. If not given, the filename will
            consist of the name of the colormap and extension '.txt'.

        """
        labels = {'RGBA':4} if self.colors.shape[1] == 4 else {'RGB': 3}
        t = Table(self.colors,labels,f'Creator: {util.execution_stamp("Colormap")}')
        t.save(self._get_file_handle(fname,'txt'))


    def save_GOM(self,fname=None):
        """
        Save as ASCII file for use in GOM Aramis.

        Parameters
        ----------
        fname : file, str, or pathlib.Path, optional
            Filename to store results. If not given, the filename will
            consist of the name of the colormap and extension '.legend'.

        """
        # ToDo: test in GOM
        GOM_str = '1 1 {name} 9 {name} '.format(name=self.name.replace(" ","_")) \
                +  '0 1 0 3 0 0 -1 9 \\ 0 0 0 255 255 255 0 0 255 ' \
                + f'30 NO_UNIT 1 1 64 64 64 255 1 0 0 0 0 0 0 3 0 {len(self.colors)}' \
                + ' '.join([f' 0 {c[0]} {c[1]} {c[2]} 255 1' for c in reversed((self.colors*255).astype(int))]) \
                + '\n'

        self._get_file_handle(fname,'legend').write(GOM_str)


    def save_gmsh(self,fname=None):
        """
        Save as ASCII file for use in gmsh.

        Parameters
        ----------
        fname : file, str, or pathlib.Path, optional
            Filename to store results. If not given, the filename will
            consist of the name of the colormap and extension '.msh'.

        """
        # ToDo: test in gmsh
        gmsh_str = 'View.ColorTable = {\n' \
                 +'\n'.join([f'{c[0]},{c[1]},{c[2]},' for c in self.colors[:,:3]*255]) \
                 +'\n}\n'
        self._get_file_handle(fname,'msh').write(gmsh_str)


    @staticmethod
    def _interpolate_msh(frac,low,high):
        """
        Interpolate in Msh color space.

        This interpolation gives a perceptually uniform colormap.

        References
        ----------
        https://www.kennethmoreland.com/color-maps/ColorMapsExpanded.pdf
        https://www.kennethmoreland.com/color-maps/diverging_map.py

        """
        def rad_diff(a,b):
            return abs(a[2]-b[2])

        def adjust_hue(msh_sat, msh_unsat):
            """If saturation of one of the two colors is much less than the other, hue of the less."""
            if msh_sat[0] >= msh_unsat[0]:
               return msh_sat[2]
            else:
                hSpin = msh_sat[1]/np.sin(msh_sat[1])*np.sqrt(msh_unsat[0]**2.0-msh_sat[0]**2)/msh_sat[0]
                if msh_sat[2] < - np.pi/3.0: hSpin *= -1.0
                return msh_sat[2] + hSpin

        lo = np.array(low)
        hi = np.array(high)

        if (lo[1] > 0.05 and hi[1] > 0.05 and rad_diff(lo,hi) > np.pi/3.0):
            M_mid = max(lo[0],hi[0],88.0)
            if frac < 0.5:
                hi = np.array([M_mid,0.0,0.0])
                frac *= 2.0
            else:
                lo = np.array([M_mid,0.0,0.0])
                frac = 2.0*frac - 1.0
        if   lo[1] < 0.05 and hi[1] > 0.05:
            lo[2] = adjust_hue(hi,lo)
        elif lo[1] > 0.05 and hi[1] < 0.05:
            hi[2] = adjust_hue(lo,hi)

        return (1.0 - frac) * lo + frac * hi


    _predefined_mpl= {'Perceptually Uniform Sequential': [
                         'viridis', 'plasma', 'inferno', 'magma', 'cividis'],
                      'Sequential': [
                         'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                         'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                         'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'],
                      'Sequential (2)': [
                         'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
                         'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
                         'hot', 'afmhot', 'gist_heat', 'copper'],
                      'Diverging': [
                         'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
                         'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'],
                      'Cyclic': ['twilight', 'twilight_shifted', 'hsv'],
                      'Qualitative': [
                         'Pastel1', 'Pastel2', 'Paired', 'Accent',
                         'Dark2', 'Set1', 'Set2', 'Set3',
                         'tab10', 'tab20', 'tab20b', 'tab20c'],
                      'Miscellaneous': [
                         'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
                         'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
                         'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar']}

    _predefined_DAMASK = {'orientation':   {'low':  [0.933334,0.878432,0.878431],
                                            'high': [0.250980,0.007843,0.000000]},
                          'strain':        {'low':  [0.941177,0.941177,0.870588],
                                            'high': [0.266667,0.266667,0.000000]},
                          'stress':        {'low':  [0.878432,0.874511,0.949019],
                                            'high': [0.000002,0.000000,0.286275]}}

    predefined = dict(**{'DAMASK':list(_predefined_DAMASK)},**_predefined_mpl)


    @staticmethod
    def _hsv2rgb(hsv):
        """
        H(ue) S(aturation) V(alue) to R(red) G(reen) B(lue).

        References
        ----------
        https://www.rapidtables.com/convert/color/hsv-to-rgb.html

        """
        sextant = np.clip(int(hsv[0]/60.),0,5)
        c = hsv[1]*hsv[2]
        x = c*(1.0 - abs((hsv[0]/60.)%2 - 1.))

        return np.array([
                         [c, x, 0],
                         [x, c, 0],
                         [0, c, x],
                         [0, x, c],
                         [x, 0, c],
                         [c, 0, x],
                        ])[sextant] + hsv[2] - c

    @staticmethod
    def _rgb2hsv(rgb):
        """
        R(ed) G(reen) B(lue) to H(ue) S(aturation) V(alue).

        References
        ----------
        https://www.rapidtables.com/convert/color/rgb-to-hsv.html

        """
        C_max = rgb.max()
        C_min = rgb.min()
        Delta = C_max - C_min

        v = C_max
        s = 0. if np.isclose(C_max,0.) else Delta/C_max
        if np.isclose(Delta,0.):
            h = 0.
        elif rgb.argmax() == 0:
            h = (rgb[1]-rgb[2])/Delta%6
        elif rgb.argmax() == 1:
            h = (rgb[2]-rgb[0])/Delta + 2.
        elif rgb.argmax() == 2:
            h = (rgb[0]-rgb[1])/Delta + 4.

        h = np.clip(h,0.,6.) * 60.

        return np.array([h,s,v])


    @staticmethod
    def _hsl2rgb(hsl):
        """
        H(ue) S(aturation) L(uminance) to R(red) G(reen) B(lue).

        References
        ----------
        https://www.rapidtables.com/convert/color/hsl-to-rgb.html

        """
        sextant = np.clip(int(hsl[0]/60.),0,5)
        c = (1.0 - abs(2.0 * hsl[2] - 1.))*hsl[1]
        x = c*(1.0 - abs((hsl[0]/60.)%2 - 1.))
        m = hsl[2] - 0.5*c

        return np.array([
                         [c+m, x+m, m],
                         [x+m, c+m, m],
                         [m, c+m, x+m],
                         [m, x+m, c+m],
                         [x+m, m, c+m],
                         [c+m, m, x+m],
                        ])[sextant]

    @staticmethod
    def _rgb2hsl(rgb):
        """
        R(ed) G(reen) B(lue) to H(ue) S(aturation) L(uminance).

        References
        ----------
        https://www.rapidtables.com/convert/color/rgb-to-hsl.html

        """
        C_max = rgb.max()
        C_min = rgb.min()
        Delta = C_max - C_min

        l = np.clip((C_max + C_min)*.5,0.,1.)                                                       # noqa
        s = 0. if np.isclose(C_max,C_min) else Delta/(1.-np.abs(2*l-1.))
        if np.isclose(Delta,0.):
            h = 0.
        elif rgb.argmax() == 0:
            h = (rgb[1]-rgb[2])/Delta%6
        elif rgb.argmax() == 1:
            h = (rgb[2]-rgb[0])/Delta + 2.
        elif rgb.argmax() == 2:
            h = (rgb[0]-rgb[1])/Delta + 4.

        h = np.clip(h,0.,6.) * 60.

        return np.array([h,s,l])


    @staticmethod
    def _xyz2rgb(xyz):
        """
        CIE Xyz to R(ed) G(reen) B(lue).

        References
        ----------
        https://www.easyrgb.com/en/math.php

        """
        rgb_lin = np.dot(np.array([
                                   [ 3.240969942,-1.537383178,-0.498610760],
                                   [-0.969243636, 1.875967502, 0.041555057],
                                   [ 0.055630080,-0.203976959, 1.056971514]
                                  ]),xyz)
        with np.errstate(invalid='ignore'):
            rgb = np.where(rgb_lin>0.0031308,rgb_lin**(1.0/2.4)*1.0555-0.0555,rgb_lin*12.92)

        return np.clip(rgb,0.,1.)

    @staticmethod
    def _rgb2xyz(rgb):
        """
        R(ed) G(reen) B(lue) to CIE Xyz.

        References
        ----------
        https://www.easyrgb.com/en/math.php

        """
        rgb_lin = np.where(rgb>0.04045,((rgb+0.0555)/1.0555)**2.4,rgb/12.92)
        return np.dot(np.array([
                                [0.412390799,0.357584339,0.180480788],
                                [0.212639006,0.715168679,0.072192315],
                                [0.019330819,0.119194780,0.950532152]
                               ]),rgb_lin)


    @staticmethod
    def _lab2xyz(lab,ref_white=None):
        """
        CIE Lab to CIE Xyz.

        References
        ----------
        http://www.brucelindbloom.com/index.html?Eqn_Lab_to_XYZ.html

        """
        f_x = (lab[0]+16.)/116. + lab[1]/500.
        f_z = (lab[0]+16.)/116. - lab[2]/200.

        return np.array([
                         f_x**3.                if f_x**3. > _eps     else (116.*f_x-16.)/_kappa,
                         ((lab[0]+16.)/116.)**3 if lab[0]>_kappa*_eps else lab[0]/_kappa,
                         f_z**3.                if f_z**3. > _eps     else (116.*f_z-16.)/_kappa
                        ])*(ref_white if ref_white is not None else _ref_white)

    @staticmethod
    def _xyz2lab(xyz,ref_white=None):
        """
        CIE Xyz to CIE Lab.

        References
        ----------
        http://www.brucelindbloom.com/index.html?Eqn_Lab_to_XYZ.html

        """
        ref_white = ref_white if ref_white is not None else _ref_white
        f = np.where(xyz/ref_white > _eps,(xyz/ref_white)**(1./3.),(_kappa*xyz/ref_white+16.)/116.)

        return np.array([
                         116.0 *  f[1] - 16.0,
                         500.0 * (f[0] - f[1]),
                         200.0 * (f[1] - f[2])
                        ])


    @staticmethod
    def _lab2msh(lab):
        """
        CIE Lab to Msh.

        References
        ----------
        https://www.kennethmoreland.com/color-maps/ColorMapsExpanded.pdf
        https://www.kennethmoreland.com/color-maps/diverging_map.py

        """
        M = np.linalg.norm(lab)
        return np.array([
                         M,
                         np.arccos(lab[0]/M)       if M>1e-8 else 0.,
                         np.arctan2(lab[2],lab[1]) if M>1e-8 else 0.,
                        ])

    @staticmethod
    def _msh2lab(msh):
        """
        Msh to CIE Lab.

        References
        ----------
        https://www.kennethmoreland.com/color-maps/ColorMapsExpanded.pdf
        https://www.kennethmoreland.com/color-maps/diverging_map.py

        """
        return np.array([
                         msh[0] * np.cos(msh[1]),
                         msh[0] * np.sin(msh[1]) * np.cos(msh[2]),
                         msh[0] * np.sin(msh[1]) * np.sin(msh[2])
                        ])

    @staticmethod
    def _lab2rgb(lab):
        return Colormap._xyz2rgb(Colormap._lab2xyz(lab))

    @staticmethod
    def _rgb2lab(rgb):
        return Colormap._xyz2lab(Colormap._rgb2xyz(rgb))

    @staticmethod
    def _msh2rgb(msh):
        return Colormap._lab2rgb(Colormap._msh2lab(msh))

    @staticmethod
    def _rgb2msh(rgb):
        return Colormap._lab2msh(Colormap._rgb2lab(rgb))

    @staticmethod
    def _hsv2msh(hsv):
        return Colormap._rgb2msh(Colormap._hsv2rgb(hsv))

    @staticmethod
    def _hsl2msh(hsl):
        return Colormap._rgb2msh(Colormap._hsl2rgb(hsl))

    @staticmethod
    def _xyz2msh(xyz):
        return Colormap._lab2msh(Colormap._xyz2lab(xyz))
