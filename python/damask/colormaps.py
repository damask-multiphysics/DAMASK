import math

import numpy as np

class Color():
  """
  Conversion of colors between different color-spaces.

  Colors should be given in the form Color('model',[vector]).
  To convert or copy color from one space to other, use the methods
  convertTo('model') or expressAs('model'), respectively.
  """

  __slots__ = [
               'model',
               'color',
               '__dict__',
              ]


# ------------------------------------------------------------------
  def __init__(self,
               model = 'RGB',
               color = np.zeros(3,'d')):

    self.__transforms__ = \
                   {'HSV':    {'index': 0, 'next': self._HSV2HSL},
                    'HSL':    {'index': 1, 'next': self._HSL2RGB,     'prev': self._HSL2HSV},
                    'RGB':    {'index': 2, 'next': self._RGB2XYZ,     'prev': self._RGB2HSL},
                    'XYZ':    {'index': 3, 'next': self._XYZ2CIELAB,  'prev': self._XYZ2RGB},
                    'CIELAB': {'index': 4, 'next': self._CIELAB2MSH,  'prev': self._CIELAB2XYZ},
                    'MSH':    {'index': 5, 'prev': self._MSH2CIELAB},
                    }

    model = model.upper()
    if model not in list(self.__transforms__.keys()): model = 'RGB'
    if model == 'RGB' and max(color) > 1.0:                                                         # are we RGB255 ?
      for i in range(3):
        color[i] /= 255.0                                                                           # rescale to RGB

    if model == 'HSL':                                                                              # are we HSL ?
      if abs(color[0]) > 1.0: color[0] /= 360.0                                                     # with angular hue?
      while color[0] >= 1.0: color[0] -= 1.0                                                        # rewind to proper range
      while color[0] < 0.0:  color[0] += 1.0                                                        # rewind to proper range

    self.model = model
    self.color = np.array(color,'d')


# ------------------------------------------------------------------
  def __repr__(self):
    """Color model and values."""
    return 'Model: %s Color: %s'%(self.model,str(self.color))


# ------------------------------------------------------------------
  def __str__(self):
    """Color model and values."""
    return self.__repr__()


# ------------------------------------------------------------------
  def convertTo(self,toModel = 'RGB'):
    toModel = toModel.upper()
    if toModel not in list(self.__transforms__.keys()): return

    sourcePos = self.__transforms__[self.model]['index']
    targetPos = self.__transforms__[toModel]['index']

    while sourcePos < targetPos:
      self.__transforms__[self.model]['next']()
      sourcePos += 1

    while sourcePos > targetPos:
      self.__transforms__[self.model]['prev']()
      sourcePos -= 1
    return self


# ------------------------------------------------------------------
  def expressAs(self,asModel = 'RGB'):
    return self.__class__(self.model,self.color).convertTo(asModel)



  def _HSV2HSL(self):
    """
    Convert H(ue) S(aturation) V(alue or brightness) to H(ue) S(aturation) L(uminance).

    All values are in the range [0,1]
    http://codeitdown.com/hsl-hsb-hsv-color/
    """
    if self.model != 'HSV': return

    converted = Color('HSL',np.array([
                                      self.color[0],
                                      1. if self.color[2] == 0.0 or (self.color[1] == 0.0 and self.color[2] == 1.0) \
                                      else self.color[1]*self.color[2]/(1.-abs(self.color[2]*(2.-self.color[1])-1.)),
                                      0.5*self.color[2]*(2.-self.color[1]),
                                     ]))

    self.model = converted.model
    self.color = converted.color


  def _HSL2HSV(self):
    """
    Convert H(ue) S(aturation) L(uminance) to H(ue) S(aturation) V(alue or brightness).

    All values are in the range [0,1]
    http://codeitdown.com/hsl-hsb-hsv-color/
    """
    if self.model != 'HSL': return

    h = self.color[0]
    b = self.color[2]+0.5*(self.color[1]*(1.-abs(2*self.color[2]-1)))
    s = 1.0 if b == 0.0 else 2.*(b-self.color[2])/b

    converted = Color('HSV',np.array([h,s,b]))

    self.model = converted.model
    self.color = converted.color


  def _HSL2RGB(self):
    """
    Convert H(ue) S(aturation) L(uminance) to R(red) G(reen) B(lue).

    All values are in the range [0,1]
    from http://en.wikipedia.org/wiki/HSL_and_HSV
    """
    if self.model != 'HSL': return

    sextant = self.color[0]*6.0
    c = (1.0 - abs(2.0 * self.color[2] - 1.0))*self.color[1]
    x = c*(1.0 - abs(sextant%2 - 1.0))
    m = self.color[2] - 0.5*c

    converted = Color('RGB',np.array([
                                      [c+m, x+m, m],
                                      [x+m, c+m, m],
                                      [m, c+m, x+m],
                                      [m, x+m, c+m],
                                      [x+m, m, c+m],
                                      [c+m, m, x+m],
                                     ][int(sextant)],'d'))
    self.model = converted.model
    self.color = converted.color


  def _RGB2HSL(self):
    """
    Convert R(ed) G(reen) B(lue) to H(ue) S(aturation) L(uminance).

    All values are in the range [0,1]
    from http://130.113.54.154/~monger/hsl-rgb.html
    """
    if self.model != 'RGB': return

    HSL = np.zeros(3,'d')
    maxcolor = self.color.max()
    mincolor = self.color.min()
    HSL[2] = (maxcolor + mincolor)/2.0
    if(mincolor == maxcolor):
      HSL[0] = 0.0
      HSL[1] = 0.0
    else:
      if (HSL[2]<0.5):
        HSL[1] = (maxcolor - mincolor)/(maxcolor + mincolor)
      else:
        HSL[1] = (maxcolor - mincolor)/(2.0 - maxcolor - mincolor)
      if (maxcolor == self.color[0]):
        HSL[0] = 0.0 + (self.color[1] - self.color[2])/(maxcolor - mincolor)
      elif (maxcolor == self.color[1]):
        HSL[0] = 2.0 + (self.color[2] - self.color[0])/(maxcolor - mincolor)
      elif (maxcolor == self.color[2]):
        HSL[0] = 4.0 + (self.color[0] - self.color[1])/(maxcolor - mincolor)
      HSL[0] = HSL[0]*60.0                                                                           # scaling to 360 might be dangerous for small values
      if (HSL[0] < 0.0):
        HSL[0] = HSL[0] + 360.0
    for i in range(2):
      HSL[i+1] = min(HSL[i+1],1.0)
      HSL[i+1] = max(HSL[i+1],0.0)

    converted = Color('HSL', HSL)
    self.model = converted.model
    self.color = converted.color



  def _RGB2XYZ(self):
    """
    Convert R(ed) G(reen) B(lue) to CIE XYZ.

    All values are in the range [0,1]
    from http://www.cs.rit.edu/~ncs/color/t_convert.html
    """
    if self.model != 'RGB': return

    XYZ     = np.zeros(3,'d')
    RGB_lin = np.zeros(3,'d')
    convert = np.array([[0.412453,0.357580,0.180423],
                        [0.212671,0.715160,0.072169],
                        [0.019334,0.119193,0.950227]])

    for i in range(3):
      if (self.color[i] > 0.04045): RGB_lin[i] = ((self.color[i]+0.0555)/1.0555)**2.4
      else:                         RGB_lin[i] =   self.color[i]        /12.92
    XYZ = np.dot(convert,RGB_lin)
    for i in range(3):

      XYZ[i] = max(XYZ[i],0.0)

    converted = Color('XYZ', XYZ)
    self.model = converted.model
    self.color = converted.color



  def _XYZ2RGB(self):
    """
    Convert  CIE XYZ to R(ed) G(reen) B(lue).

    All values are in the range [0,1]
    from http://www.cs.rit.edu/~ncs/color/t_convert.html
    """
    if self.model != 'XYZ':
      return

    convert = np.array([[ 3.240479,-1.537150,-0.498535],
                        [-0.969256, 1.875992, 0.041556],
                        [ 0.055648,-0.204043, 1.057311]])
    RGB_lin = np.dot(convert,self.color)
    RGB     = np.zeros(3,'d')

    for i in range(3):
      if (RGB_lin[i] > 0.0031308): RGB[i] = ((RGB_lin[i])**(1.0/2.4))*1.0555-0.0555
      else:                        RGB[i] =   RGB_lin[i]             *12.92
    for i in range(3):
      RGB[i] = min(RGB[i],1.0)
      RGB[i] = max(RGB[i],0.0)

    maxVal = max(RGB)                                                                               # clipping colors according to the display gamut
    if (maxVal > 1.0): RGB /= maxVal

    converted = Color('RGB', RGB)
    self.model = converted.model
    self.color = converted.color



  def _CIELAB2XYZ(self):
    """
    Convert  CIE Lab to CIE XYZ.

    All values are in the range [0,1]
    from http://www.easyrgb.com/index.php?X=MATH&H=07#text7
    """
    if self.model != 'CIELAB': return

    ref_white = np.array([.95047, 1.00000, 1.08883],'d')                                            # Observer = 2, Illuminant = D65
    XYZ       = np.zeros(3,'d')

    XYZ[1] = (self.color[0] + 16.0 ) / 116.0
    XYZ[0] = XYZ[1]   + self.color[1]/ 500.0
    XYZ[2] = XYZ[1]   - self.color[2]/ 200.0

    for i in range(len(XYZ)):
      if (XYZ[i] > 6./29. ): XYZ[i] = XYZ[i]**3.
      else:                  XYZ[i] = 108./841. * (XYZ[i] - 4./29.)

    converted = Color('XYZ', XYZ*ref_white)
    self.model = converted.model
    self.color = converted.color

  def _XYZ2CIELAB(self):
    """
    Convert CIE XYZ to CIE Lab.

    All values are in the range [0,1]
    from http://en.wikipedia.org/wiki/Lab_color_space,
         http://www.cs.rit.edu/~ncs/color/t_convert.html
    """
    if self.model != 'XYZ': return

    ref_white = np.array([.95047, 1.00000, 1.08883],'d')                                            # Observer = 2, Illuminant = D65
    XYZ = self.color/ref_white

    for i in range(len(XYZ)):
      if (XYZ[i] > 216./24389 ): XYZ[i] = XYZ[i]**(1.0/3.0)
      else:                      XYZ[i] = (841./108. * XYZ[i]) + 16.0/116.0

    converted = Color('CIELAB', np.array([ 116.0 *  XYZ[1] - 16.0,
                                           500.0 * (XYZ[0] - XYZ[1]),
                                           200.0 * (XYZ[1] - XYZ[2]) ]))
    self.model = converted.model
    self.color = converted.color


  def _CIELAB2MSH(self):
    """
    Convert CIE Lab to Msh colorspace.

    from http://www.cs.unm.edu/~kmorel/documents/ColorMaps/DivergingColorMapWorkshop.xls
    """
    if self.model != 'CIELAB': return

    Msh = np.zeros(3,'d')
    Msh[0] = math.sqrt(np.dot(self.color,self.color))
    if (Msh[0] > 0.001):
      Msh[1] = math.acos(self.color[0]/Msh[0])
      if (self.color[1] != 0.0):
        Msh[2] = math.atan2(self.color[2],self.color[1])

    converted = Color('MSH', Msh)
    self.model = converted.model
    self.color = converted.color


  def _MSH2CIELAB(self):
    """
    Convert Msh colorspace to CIE Lab.

    with s,h in radians
    from http://www.cs.unm.edu/~kmorel/documents/ColorMaps/DivergingColorMapWorkshop.xls
    """
    if self.model != 'MSH': return

    Lab = np.zeros(3,'d')
    Lab[0] = self.color[0] * math.cos(self.color[1])
    Lab[1] = self.color[0] * math.sin(self.color[1]) * math.cos(self.color[2])
    Lab[2] = self.color[0] * math.sin(self.color[1]) * math.sin(self.color[2])

    converted = Color('CIELAB', Lab)
    self.model = converted.model
    self.color = converted.color


class Colormap():
  """Perceptually uniform diverging or sequential colormaps."""

  __slots__ = [
               'left',
               'right',
               'interpolate',
              ]
  __predefined__ = {
                     'gray':          {'left':  Color('HSL',[0,1,1]),
                                       'right': Color('HSL',[0,0,0.15]),
                                       'interpolate': 'perceptualuniform'},
                     'grey':          {'left':  Color('HSL',[0,1,1]),
                                       'right': Color('HSL',[0,0,0.15]),
                                       'interpolate': 'perceptualuniform'},
                     'red':           {'left':  Color('HSL',[0,1,0.14]),
                                       'right': Color('HSL',[0,0.35,0.91]),
                                       'interpolate': 'perceptualuniform'},
                     'green':         {'left':  Color('HSL',[0.33333,1,0.14]),
                                       'right': Color('HSL',[0.33333,0.35,0.91]),
                                       'interpolate': 'perceptualuniform'},
                     'blue':          {'left':  Color('HSL',[0.66,1,0.14]),
                                       'right': Color('HSL',[0.66,0.35,0.91]),
                                       'interpolate': 'perceptualuniform'},
                     'seaweed':       {'left':  Color('HSL',[0.78,1.0,0.1]),
                                       'right': Color('HSL',[0.40000,0.1,0.9]),
                                       'interpolate': 'perceptualuniform'},
                     'bluebrown':     {'left':  Color('HSL',[0.65,0.53,0.49]),
                                       'right': Color('HSL',[0.11,0.75,0.38]),
                                       'interpolate': 'perceptualuniform'},
                     'redgreen':      {'left':  Color('HSL',[0.97,0.96,0.36]),
                                       'right': Color('HSL',[0.33333,1.0,0.14]),
                                       'interpolate': 'perceptualuniform'},
                     'bluered':       {'left':  Color('HSL',[0.65,0.53,0.49]),
                                       'right': Color('HSL',[0.97,0.96,0.36]),
                                       'interpolate': 'perceptualuniform'},
                     'blueredrainbow':{'left':  Color('HSL',[2.0/3.0,1,0.5]),
                                       'right': Color('HSL',[0,1,0.5]),
                                       'interpolate': 'linear'           },
                     'orientation':   {'left':  Color('RGB',[0.933334,0.878432,0.878431]),
                                       'right': Color('RGB',[0.250980,0.007843,0.000000]),
                                       'interpolate': 'perceptualuniform'},
                     'strain':        {'left':  Color('RGB',[0.941177,0.941177,0.870588]),
                                       'right': Color('RGB',[0.266667,0.266667,0.000000]),
                                       'interpolate': 'perceptualuniform'},
                     'stress':        {'left':  Color('RGB',[0.878432,0.874511,0.949019]),
                                       'right': Color('RGB',[0.000002,0.000000,0.286275]),
                                       'interpolate': 'perceptualuniform'},
                    }


# ------------------------------------------------------------------
  def __init__(self,
               left  = Color('RGB',[1,1,1]),
               right = Color('RGB',[0,0,0]),
               interpolate = 'perceptualuniform',
               predefined = None
               ):

    if predefined is not None:
      left = self.__predefined__[predefined.lower()]['left']
      right= self.__predefined__[predefined.lower()]['right']
      interpolate = self.__predefined__[predefined.lower()]['interpolate']

    if left.__class__.__name__ != 'Color':
      left = Color()
    if right.__class__.__name__ != 'Color':
      right = Color()

    self.left  = left
    self.right = right
    self.interpolate = interpolate


# ------------------------------------------------------------------
  def __repr__(self):
    """Left and right value of colormap."""
    return 'Left: %s Right: %s'%(self.left,self.right)


# ------------------------------------------------------------------
  def invert(self):
    (self.left, self.right) = (self.right, self.left)
    return self


# ------------------------------------------------------------------
  def color(self,fraction = 0.5):

    def interpolate_Msh(lo, hi, frac):

      def rad_diff(a,b):
        return abs(a[2]-b[2])
# if saturation of one of the two colors is too less than the other, hue of the less
      def adjust_hue(Msh_sat, Msh_unsat):
        if Msh_sat[0] >= Msh_unsat[0]:
          return Msh_sat[2]
        else:
          hSpin = Msh_sat[1]/math.sin(Msh_sat[1])*math.sqrt(Msh_unsat[0]**2.0-Msh_sat[0]**2)/Msh_sat[0]
          if Msh_sat[2] < - math.pi/3.0: hSpin *= -1.0
          return Msh_sat[2] + hSpin

      Msh1 = np.array(lo[:])
      Msh2 = np.array(hi[:])

      if (Msh1[1] > 0.05 and Msh2[1] > 0.05 and rad_diff(Msh1,Msh2) > math.pi/3.0):
        M_mid = max(Msh1[0],Msh2[0],88.0)
        if frac < 0.5:
          Msh2 = np.array([M_mid,0.0,0.0],'d')
          frac *= 2.0
        else:
          Msh1 = np.array([M_mid,0.0,0.0],'d')
          frac = 2.0*frac - 1.0
      if   Msh1[1] < 0.05 and Msh2[1] > 0.05: Msh1[2] = adjust_hue(Msh2,Msh1)
      elif Msh1[1] > 0.05 and Msh2[1] < 0.05: Msh2[2] = adjust_hue(Msh1,Msh2)
      Msh = (1.0 - frac) * Msh1 + frac * Msh2

      return Color('MSH',Msh)

    def interpolate_linear(lo, hi, frac):
      """Linear interpolation between lo and hi color at given fraction; output in model of lo color."""
      interpolation = (1.0 - frac) * np.array(lo.color[:]) \
                           + frac  * np.array(hi.expressAs(lo.model).color[:])

      return Color(lo.model,interpolation)

    if self.interpolate == 'perceptualuniform':
      return interpolate_Msh(self.left.expressAs('MSH').color,
                             self.right.expressAs('MSH').color,fraction)
    elif self.interpolate == 'linear':
      return interpolate_linear(self.left,
                                self.right,fraction)
    else:
      raise NameError('unknown color interpolation method')

# ------------------------------------------------------------------
  def export(self,name = 'uniformPerceptualColorMap',\
                  format = 'paraview',\
                  steps = 2,\
                  crop = [-1.0,1.0],
                  model = 'RGB'):
    """
    [RGB] colormap for use in paraview or gmsh, or as raw string, or array.

    Arguments: name, format, steps, crop.
    Format is one of (paraview, gmsh, raw, list).
    Crop selects a (sub)range in [-1.0,1.0].
    Generates sequential map if one limiting color is either white or black,
    diverging map otherwise.
    """
    format = format.lower()                                                                         # consistent comparison basis
    frac = 0.5*(np.array(crop) + 1.0)                                                               # rescale crop range to fractions
    colors = [self.color(float(i)/(steps-1)*(frac[1]-frac[0])+frac[0]).expressAs(model).color for i in range(steps)]
    if   format == 'paraview':
      colormap = ['[\n {{\n  "ColorSpace": "RGB", "Name": "{}", "DefaultMap": true,\n  "RGBPoints" : ['.format(name)] \
               + ['    {:4d},{:8.6f},{:8.6f},{:8.6f},'.format(i,color[0],color[1],color[2],)
                                                                        for i,color in enumerate(colors[:-1])]\
               + ['    {:4d},{:8.6f},{:8.6f},{:8.6f} '.format(len(colors),colors[-1][0],colors[-1][1],colors[-1][2],)]\
               + ['   ]\n }\n]']
    elif format == 'gmsh':
      colormap = ['View.ColorTable = {'] \
               + [',\n'.join(['{%s}'%(','.join([str(x*255.0) for x in color])) for color in colors])] \
               + ['}']

    elif format == 'gom':
      colormap = ['1 1 ' + str(name)
                 + ' 9 ' + str(name)
                 + ' 0 1 0 3 0 0 -1 9 \ 0 0 0 255 255 255 0 0 255 '
                 + '30 NO_UNIT 1 1 64 64 64 255 1 0 0 0 0 0 0 3 0 ' + str(len(colors))
                 + ' '.join([' 0 %s 255 1'%(' '.join([str(int(x*255.0)) for x in color])) for color in reversed(colors)])]

    elif format == 'raw':
      colormap = ['\t'.join(map(str,color)) for color in colors]

    elif format == 'list':
      colormap = colors

    else:
      raise NameError('unknown color export format')

    return '\n'.join(colormap) + '\n' if type(colormap[0]) is str else colormap
