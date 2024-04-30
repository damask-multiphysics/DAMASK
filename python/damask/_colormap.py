import os
import json
import functools
import colorsys
from typing import Optional, Union
from itertools import chain

import numpy as np
import scipy.interpolate as interp
import matplotlib as mpl
if os.name == 'posix' and 'DISPLAY' not in os.environ:
    mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from PIL import Image

from ._typehints import FloatSequence, FileHandle
from . import util
from . import Table

_EPS   = 216./24389.
_KAPPA = 24389./27.
_REF_WHITE = np.array([.95047, 1.00000, 1.08883])                                                   # Observer = 2, Illuminant = D65

# ToDo (if needed)
# - support alpha channel (paraview/ASCII/input)
# - support NaN color (paraview)

class Colormap(mpl.colors.ListedColormap):
    """
    Enhance matplotlib colormap functionality for use within DAMASK.

    Colors are internally stored as R(ed) G(green) B(lue) values.
    A colormap can be used in matplotlib, seaborn, etc., or can be
    exported to file for external use.

    References
    ----------
    K. Moreland, Proceedings of the 5th International Symposium on Advances in Visual Computing, 2009
    https://doi.org/10.1007/978-3-642-10520-3_9

    P. Eisenlohr et al., International Journal of Plasticity 46:37–53, 2013
    https://doi.org/10.1016/j.ijplas.2012.09.012

    Matplotlib colormaps overview
    https://matplotlib.org/stable/tutorials/colors/colormaps.html

    """

    def __init__(self,
                 colors: np.ndarray, name: str):
        """
        New colormap.

        Parameters
        ----------
        colors : numpy.ndarray, shape (:,3) or (:,4)
            Color specifications as RGB(A) values.
        name : str
            String to identify the colormap.

        """
        super().__init__(colors,name)
        self.colors: np.ndarray = np.asarray(colors)


    def __eq__(self,
               other: object) -> bool:
        """
        Return self==other.

        Test equality of other.

        """
        if not isinstance(other, Colormap):
            return NotImplemented
        return np.array_equal(self.colors,other.colors)


    def __add__(self,
                other: 'Colormap') -> 'Colormap':
        """
        Return self+other.

        Concatenate.

        """
        return Colormap(np.vstack((self.colors,other.colors)),
                        f'{self.name}+{other.name}')


    def __iadd__(self,
                 other: 'Colormap') -> 'Colormap':
        """
        Return self+=other.

        Concatenate (in-place).

        """
        return self.__add__(other)


    def __mul__(self,
                factor: int) -> 'Colormap':
        """
        Return self*other.

        Repeat.

        """
        return Colormap(np.vstack([self.colors]*factor),f'{self.name}*{factor}')

    def __imul__(self,
                 factor: int) -> 'Colormap':
        """
        Return self*=other.

        Repeat (in-place).

        """
        return self.__mul__(factor)


    def __invert__(self) -> 'Colormap':
        """
        Return ~self.

        Reverse.

        """
        return self.reversed()


    def __repr__(self) -> str:
        """
        Return repr(self).

        Show as matplotlib figure.

        """
        fig = plt.figure(self.name,figsize=(5,.5))
        ax1 = fig.add_axes((0, 0, 1, 1))
        ax1.set_axis_off()
        ax1.imshow(np.linspace(0,1,self.N).reshape(1,-1),
                   aspect='auto', cmap=self, interpolation='nearest')
        plt.show(block=False)
        return f'Colormap: {self.name}'


    @staticmethod
    def from_range(low: FloatSequence,
                   high: FloatSequence,
                   name: str = 'DAMASK colormap',
                   N: int = 256,
                   model: str = 'rgb') -> 'Colormap':
        """
        Create a perceptually uniform colormap between given (inclusive) bounds.

        Parameters
        ----------
        low : sequence of float, len (3)
            Color definition for minimum value.
        high : sequence of float, len (3)
            Color definition for maximum value.
        name : str, optional
            Name of the colormap. Defaults to 'DAMASK colormap'.
        N : int, optional
            Number of color quantization levels. Defaults to 256.
        model : {'rgb', 'hsv', 'hsl', 'xyz', 'lab', 'msh'}
            Color model used for input color definitions. Defaults to 'rgb'.
            The available color models are:
            - 'rgb': Red Green Blue.
            - 'hsv': Hue Saturation Value.
            - 'hsl': Hue Saturation Luminance.
            - 'xyz': CIE Xyz.
            - 'lab': CIE Lab.
            - 'msh': Msh (for perceptually uniform interpolation).

        Returns
        -------
        new : damask.Colormap
            Colormap spanning given bounds.

        Examples
        --------
        >>> import damask
        >>> damask.Colormap.from_range((0,0,1),(0,0,0),'blue_to_black')

        """
        toMsh = dict(
            rgb=Colormap._rgb2msh,
            hsv=Colormap._hsv2msh,
            hsl=Colormap._hsl2msh,
            xyz=Colormap._xyz2msh,
            lab=Colormap._lab2msh,
            msh=lambda x:x,
        )

        if model.lower() not in toMsh:
            raise ValueError(f'invalid color model "{model}"')

        low_high = np.vstack((low,high)).astype(float)
        out_of_bounds = np.bool_(False)

        if   model.lower() == 'rgb':
            out_of_bounds = np.any(low_high<0) or np.any(low_high>1)
        elif model.lower() == 'hsv':
            out_of_bounds = np.any(low_high<0) or np.any(low_high>[360,1,1])
        elif model.lower() == 'hsl':
            out_of_bounds = np.any(low_high<0) or np.any(low_high>[360,1,1])
        elif model.lower() == 'lab':
            out_of_bounds = np.any(low_high[:,0]<0)

        if out_of_bounds:
            raise ValueError(f'{model.upper()} colors {low_high[0]} | {low_high[1]} are out of bounds')

        low_,high_ = map(toMsh[model.lower()],low_high)
        msh = map(functools.partial(Colormap._interpolate_msh,low=low_,high=high_),np.linspace(0,1,N))
        rgb = np.array(list(map(Colormap._msh2rgb,msh)))

        return Colormap(rgb,name=name)


    @staticmethod
    def from_predefined(name: str,
                        N: int = 256) -> 'Colormap':
        """
        Select from a set of predefined colormaps.

        Predefined colormaps (Colormap.predefined) include
        native matplotlib colormaps and common DAMASK colormaps.

        Parameters
        ----------
        name : str
            Name of the colormap.
        N : int, optional
            Number of color quantization levels. Defaults to 256.
            This parameter is not used for matplotlib colormaps
            that are of type `ListedColormap`.

        Returns
        -------
        new : damask.Colormap
            Predefined colormap.

        Examples
        --------
        >>> import damask
        >>> damask.Colormap.from_predefined('strain')

        """
        try:
            # matplotlib presets
            colormap = cm.__dict__[name]
            return Colormap(np.array(list(map(colormap,np.linspace(0,1,N)))
                                     if isinstance(colormap,mpl.colors.LinearSegmentedColormap) else
                                     colormap.colors),
                            name=name)
        except KeyError:
            # DAMASK presets
            definition = Colormap._predefined_DAMASK[name]
            return Colormap.from_range(definition['low'],definition['high'],name,N)


    def at(self,
           fraction : Union[float,FloatSequence]) -> np.ndarray:
        """
        Interpolate color at fraction.

        Parameters
        ----------
        fraction : (sequence of) float
            Fractional coordinate(s) to evaluate Colormap at.

        Returns
        -------
        color : numpy.ndarray, shape(...,4)
            RGBA values of interpolated color(s).

        Examples
        --------
        >>> import damask
        >>> cmap = damask.Colormap.from_predefined('gray')
        >>> cmap.at(0.5)
        array([0.5, 0.5, 0.5, 1. ])
        >>> 'rgb({},{},{})'.format(*cmap.at(0.5))
        'rgb(0.5,0.5,0.5)'

        """
        return interp.interp1d(np.linspace(0,1,self.N),
                               self.colors,
                               axis=0,
                               assume_sorted=True)(fraction)


    def shade(self,
              field: np.ndarray,
              bounds: Optional[FloatSequence] = None,
              gap: Optional[float] = None) -> Image.Image:
        """
        Generate PIL image of 2D field using colormap.

        Parameters
        ----------
        field : numpy.ndarray, shape (:,:)
            Data to be shaded.
        bounds : sequence of float, len (2), optional
            Value range (left,right) spanned by colormap.
        gap : field.dtype, optional
            Transparent value. NaN will always be rendered transparent.
            Defaults to None.

        Returns
        -------
        PIL.Image
            RGBA image of shaded data.

        """
        mask = np.logical_not(np.isnan(field) if gap is None else
                              np.logical_or(np.isnan(field), field == gap))                         # mask NaN (and gap if present)

        l,r = (field[mask].min(),field[mask].max()) if bounds is None else \
              (bounds[0],bounds[1])

        if abs(delta := r-l) * 1e8 <= (avg := 0.5*abs(r+l)):                                        # delta is similar to numerical noise
            l,r = (l-0.5*avg*np.sign(delta),r+0.5*avg*np.sign(delta))                               # extend range to have actual data centered within

        field_ = np.nan_to_num(field, nan=(l+r)/2, posinf=r, neginf=l)

        return Image.fromarray(
            (np.dstack((
                        self.colors[np.round(np.clip((field_-l)/(r-l),0.0,1.0)*(self.N-1)).astype(np.uint16),:3],
                        mask.astype(float)
                       )
                      )*255
            ).astype(np.uint8),
            mode='RGBA')


    def reversed(self,
                 name: Optional[str] = None) -> 'Colormap':
        """
        Reverse.

        Parameters
        ----------
        name : str, optional
            Name of the reversed colormap.
            Defaults to parent colormap name + '_r'.

        Returns
        -------
        damask.Colormap
            Reversed colormap.

        Examples
        --------
        >>> import damask
        >>> damask.Colormap.from_predefined('stress').reversed()
        Colormap: stress_r

        """
        rev = super().reversed(name)
        return Colormap(np.array(rev.colors),rev.name[:-4] if rev.name.endswith('_r_r') else rev.name)


    def save_paraview(self,
                      fname: Optional[FileHandle] = None):
        """
        Save as JSON file for use in Paraview.

        Parameters
        ----------
        fname : file, str, or pathlib.Path, optional
            File to store results. Defaults to colormap name + '.json'.

        """
        out = [{
                'Creator':util.execution_stamp('Colormap'),
                'ColorSpace':'RGB',
                'Name':self.name,
                'DefaultMap':True,
                'RGBPoints':list(chain.from_iterable([(i,*c) for i,c in enumerate(self.colors.round(6))]))
               }]

        with util.open_text(self.name.replace(' ','_')+'.json' if fname is None else fname, 'w') as fhandle:
            json.dump(out,fhandle,indent=4)
            fhandle.write('\n')


    def save_ASCII(self,
                   fname: Optional[FileHandle] = None):
        """
        Save as ASCII file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path, optional
            File to store results. Defaults to colormap name + '.txt'.

        """
        labels = {'RGBA':4} if self.colors.shape[1] == 4 else {'RGB': 3}
        t = Table(labels,self.colors,[f'Creator: {util.execution_stamp("Colormap")}'])

        with util.open_text(self.name.replace(' ','_')+'.txt' if fname is None else fname, 'w') as fhandle:
            t.save(fhandle)


    def save_GOM(self, fname: Optional[FileHandle] = None):
        """
        Save as ASCII file for use in GOM Aramis.

        Parameters
        ----------
        fname : file, str, or pathlib.Path, optional
            File to store results. Defaults to colormap name + '.legend'.

        """
        # ToDo: test in GOM
        GOM_str = '1 1 {name} 9 {name} '.format(name=self.name.replace(" ","_")) \
                +  '0 1 0 3 0 0 -1 9 \\ 0 0 0 255 255 255 0 0 255 ' \
                + f'30 NO_UNIT 1 1 64 64 64 255 1 0 0 0 0 0 0 3 0 {self.N}' \
                + ' '.join([f' 0 {c[0]} {c[1]} {c[2]} 255 1' for c in reversed((self.colors*255).astype(np.int64))]) \
                + '\n'

        with util.open_text(self.name.replace(' ','_')+'.legend' if fname is None else fname, 'w') as fhandle:
            fhandle.write(GOM_str)


    def save_gmsh(self,
                  fname: Optional[FileHandle] = None):
        """
        Save as ASCII file for use in gmsh.

        Parameters
        ----------
        fname : file, str, or pathlib.Path, optional
            File to store results. Defaults to colormap name + '.msh'.

        """
        # ToDo: test in gmsh
        gmsh_str = 'View.ColorTable = {\n' \
                 +'\n'.join([f'{c[0]},{c[1]},{c[2]},' for c in self.colors[:,:3]*255]) \
                 +'\n}\n'

        with util.open_text(self.name.replace(' ','_')+'.msh' if fname is None else fname, 'w') as fhandle:
            fhandle.write(gmsh_str)


    @staticmethod
    def _interpolate_msh(frac: float,
                         low: np.ndarray,
                         high: np.ndarray) -> np.ndarray:
        """
        Interpolate in Msh color space.

        This interpolation gives a perceptually uniform colormap.

        References
        ----------
        | https://www.kennethmoreland.com/color-maps/ColorMapsExpanded.pdf
        | https://www.kennethmoreland.com/color-maps/diverging_map.py

        """
        def rad_diff(a,b):
            return abs(a[2]-b[2])

        def adjust_hue(msh_sat, msh_unsat):
            """If saturation of one of the two colors is much less than the other, hue of the less."""
            if msh_sat[0] >= msh_unsat[0]:
                return msh_sat[2]

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
        if   lo[1] < 0.05 < hi[1]:
            lo[2] = adjust_hue(hi,lo)
        elif hi[1] < 0.05 < lo[1]:
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

    _predefined_DAMASK = {'orientation':   {'low':  [0.933334,0.878432,0.878431],                   # noqa
                                            'high': [0.250980,0.007843,0.000000]},
                          'strain':        {'low':  [0.941177,0.941177,0.870588],
                                            'high': [0.266667,0.266667,0.000000]},
                          'stress':        {'low':  [0.878432,0.874511,0.949019],
                                            'high': [0.000002,0.000000,0.286275]}}

    predefined = dict(**{'DAMASK':list(_predefined_DAMASK)},**_predefined_mpl)


    @staticmethod
    def _hsv2rgb(hsv: np.ndarray) -> np.ndarray:
        """
        Hue Saturation Value to Red Green Blue.

        Parameters
        ----------
        hsv : numpy.ndarray, shape (3)
            HSV values.

        Returns
        -------
        rgb : numpy.ndarray, shape (3)
            RGB values.

        """
        return np.array(colorsys.hsv_to_rgb(hsv[0]/360.,hsv[1],hsv[2]))

    @staticmethod
    def _rgb2hsv(rgb: np.ndarray) -> np.ndarray:
        """
        Red Green Blue to Hue Saturation Value.

        Parameters
        ----------
        rgb : numpy.ndarray, shape (3)
            RGB values.

        Returns
        -------
        hsv : numpy.ndarray, shape (3)
            HSV values.

        """
        h,s,v = colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2])
        return np.array([h*360,s,v])


    @staticmethod
    def _hsl2rgb(hsl: np.ndarray) -> np.ndarray:
        """
        Hue Saturation Luminance to Red Green Blue.

        Parameters
        ----------
        hsl : numpy.ndarray, shape (3)
            HSL values.

        Returns
        -------
        rgb : numpy.ndarray, shape (3)
            RGB values.

        """
        return np.array(colorsys.hls_to_rgb(hsl[0]/360.,hsl[2],hsl[1]))

    @staticmethod
    def _rgb2hsl(rgb: np.ndarray) -> np.ndarray:
        """
        Red Green Blue to Hue Saturation Luminance.

        Parameters
        ----------
        rgb : numpy.ndarray, shape (3)
            RGB values.

        Returns
        -------
        hsl : numpy.ndarray, shape (3)
            HSL values.

        """
        h,l,s = colorsys.rgb_to_hls(rgb[0],rgb[1],rgb[2])
        return np.array([h*360,s,l])


    @staticmethod
    def _xyz2rgb(xyz: np.ndarray) -> np.ndarray:
        """
        CIE Xyz to Red Green Blue.

        Parameters
        ----------
        xyz : numpy.ndarray, shape (3)
            CIE Xyz values.

        Returns
        -------
        rgb : numpy.ndarray, shape (3)
            RGB values.

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
    def _rgb2xyz(rgb: np.ndarray) -> np.ndarray:
        """
        Red Green Blue to CIE Xyz.

        Parameters
        ----------
        rgb : numpy.ndarray, shape (3)
            RGB values.

        Returns
        -------
        xyz : numpy.ndarray, shape (3)
            CIE Xyz values.

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
    def _lab2xyz(lab: np.ndarray,
                 ref_white: np.ndarray = _REF_WHITE) -> np.ndarray:
        """
        CIE Lab to CIE Xyz.

        Parameters
        ----------
        lab : numpy.ndarray, shape (3)
            CIE lab values.
        ref_white : numpy.ndarray, shape (3)
            Reference white, default value is the standard 2° observer for D65.

        Returns
        -------
        xyz : numpy.ndarray, shape (3)
            CIE Xyz values.

        References
        ----------
        http://www.brucelindbloom.com/index.html?Eqn_Lab_to_XYZ.html

        """
        f_x = (lab[0]+16.)/116. + lab[1]/500.
        f_z = (lab[0]+16.)/116. - lab[2]/200.

        return np.array([
                         f_x**3.                if f_x**3. > _EPS     else (116.*f_x-16.)/_KAPPA,
                         ((lab[0]+16.)/116.)**3 if lab[0]>_KAPPA*_EPS else lab[0]/_KAPPA,
                         f_z**3.                if f_z**3. > _EPS     else (116.*f_z-16.)/_KAPPA
                        ])*ref_white

    @staticmethod
    def _xyz2lab(xyz: np.ndarray,
                 ref_white: np.ndarray = _REF_WHITE) -> np.ndarray:
        """
        CIE Xyz to CIE Lab.

        Parameters
        ----------
        xyz : numpy.ndarray, shape (3)
            CIE Xyz values.
        ref_white : numpy.ndarray, shape (3)
            Reference white, default value is the standard 2° observer for D65.

        Returns
        -------
        lab : numpy.ndarray, shape (3)
            CIE lab values.

        References
        ----------
        http://www.brucelindbloom.com/index.html?Eqn_Lab_to_XYZ.html

        """
        f = np.where(xyz/ref_white > _EPS,(xyz/ref_white)**(1./3.),(_KAPPA*xyz/ref_white+16.)/116.)

        return np.array([
                         116.0 *  f[1] - 16.0,
                         500.0 * (f[0] - f[1]),
                         200.0 * (f[1] - f[2])
                        ])


    @staticmethod
    def _lab2msh(lab: np.ndarray) -> np.ndarray:
        """
        CIE Lab to Msh.

        Parameters
        ----------
        lab : numpy.ndarray, shape (3)
            CIE lab values.

        Returns
        -------
        msh : numpy.ndarray, shape (3)
            Msh values.

        References
        ----------
        | https://www.kennethmoreland.com/color-maps/ColorMapsExpanded.pdf
        | https://www.kennethmoreland.com/color-maps/diverging_map.py

        """
        M = np.linalg.norm(lab)
        return np.array([
                         M,
                         np.arccos(lab[0]/M)       if M>1e-8 else 0.,
                         np.arctan2(lab[2],lab[1]) if M>1e-8 else 0.,
                        ])

    @staticmethod
    def _msh2lab(msh: np.ndarray) -> np.ndarray:
        """
        Msh to CIE Lab.

        Parameters
        ----------
        msh : numpy.ndarray, shape (3)
            Msh values.

        Returns
        -------
        lab : numpy.ndarray, shape (3)
            CIE lab values.

        References
        ----------
        | https://www.kennethmoreland.com/color-maps/ColorMapsExpanded.pdf
        | https://www.kennethmoreland.com/color-maps/diverging_map.py

        """
        return np.array([
                         msh[0] * np.cos(msh[1]),
                         msh[0] * np.sin(msh[1]) * np.cos(msh[2]),
                         msh[0] * np.sin(msh[1]) * np.sin(msh[2])
                        ])

    @staticmethod
    def _lab2rgb(lab: np.ndarray) -> np.ndarray:
        return Colormap._xyz2rgb(Colormap._lab2xyz(lab))

    @staticmethod
    def _rgb2lab(rgb: np.ndarray) -> np.ndarray:
        return Colormap._xyz2lab(Colormap._rgb2xyz(rgb))

    @staticmethod
    def _msh2rgb(msh: np.ndarray) -> np.ndarray:
        return Colormap._lab2rgb(Colormap._msh2lab(msh))

    @staticmethod
    def _rgb2msh(rgb: np.ndarray) -> np.ndarray:
        return Colormap._lab2msh(Colormap._rgb2lab(rgb))

    @staticmethod
    def _hsv2msh(hsv: np.ndarray) -> np.ndarray:
        return Colormap._rgb2msh(Colormap._hsv2rgb(hsv))

    @staticmethod
    def _hsl2msh(hsl: np.ndarray) -> np.ndarray:
        return Colormap._rgb2msh(Colormap._hsl2rgb(hsl))

    @staticmethod
    def _xyz2msh(xyz: np.ndarray) -> np.ndarray:
        return Colormap._lab2msh(Colormap._xyz2lab(xyz))
