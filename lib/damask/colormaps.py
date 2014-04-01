# -*- coding: UTF-8 no BOM -*-

### --- COLOR CLASS --------------------------------------------------

class Color():
  '''
    Conversion of colors between different color-spaces. Colors should be given in the form 
    Color('model',[vector]).To convert and copy color from one space to other, use the methods 
    convertTo('model') and expressAs('model')spectively
  '''
  import numpy
  
  __slots__ = [
               'model',
               'color',
              ]

  
  # ------------------------------------------------------------------
  def __init__(self,
               model = 'RGB',
               color = numpy.zeros(3,'d')):
    import numpy
    
    self.__transforms__ = \
                   {'HSL':    {'index': 0, 'next': self._HSL2RGB},
                    'RGB':    {'index': 1, 'next': self._RGB2XYZ,     'prev': self._RGB2HSL},
                    'XYZ':    {'index': 2, 'next': self._XYZ2CIELAB,  'prev': self._XYZ2RGB},
                    'CIELAB': {'index': 3, 'next': self._CIELAB2MSH,  'prev': self._CIELAB2XYZ},
                    'MSH':    {'index': 4, 'prev': self._MSH2CIELAB},
                    }

    model = model.upper()
    if model not in self.__transforms__.keys(): model = 'RGB'
    if model == 'RGB' and max(color) > 1.0:                                                                                                              # are we RGB255 ?
      for i in range(3):
        color[i] /= 255.0                                                                                                                                 # rescale to RGB

    if model == 'HSL':                                                                                                                                   # are we HSL ?
      if abs(color[0]) > 1.0: color[0] /= 360.0                                                                                                          # with angular hue?
      while color[0] >= 1.0: color[0] -= 1.0                                                                                                             # rewind to proper range
      while color[0] < 0.0:  color[0] += 1.0                                                                                                             # rewind to proper range

    self.model = model
    self.color = numpy.array(color,'d')
  
  
  # ------------------------------------------------------------------
  def __repr__(self):
    return 'Model: %s Color: %s'%(self.model,str(self.color))    
  
  
  # ------------------------------------------------------------------
  def __str__(self):
    return self.__repr__()
  
  
  # ------------------------------------------------------------------
  def convertTo(self,toModel = 'RGB'):
    toModel = toModel.upper()
    if toModel not in self.__transforms__.keys(): return 

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
  
  
  # ------------------------------------------------------------------   
  # convert H(ue) S(aturation) L(uminance) to R(red) G(reen) B(lue) 
  # with S,L,H,R,G,B running from 0 to 1
  # from http://en.wikipedia.org/wiki/HSL_and_HSV
  def _HSL2RGB(self):
    import numpy
    if self.model != 'HSL': return
    
    sextant = self.color[0]*6.0
    c = (1.0 - abs(2.0 * self.color[2] - 1.0))*self.color[1]
    x = c*(1.0 - abs(sextant%2 - 1.0))
    m = self.color[2] - 0.5*c

    converted = Color('RGB',numpy.array([
                                    [c+m, x+m, m],
                                    [x+m, c+m, m],
                                    [m, c+m, x+m],
                                    [m, x+m, c+m],
                                    [x+m, m, c+m],
                                    [c+m, m, x+m],
                                   ][int(sextant)],'d'))
    self.model = converted.model
    self.color = converted.color
  
  
  # ------------------------------------------------------------------
  # convert R(ed) G(reen) B(lue) to H(ue) S(aturation) L(uminance)
  # with S,L,H,R,G,B running from 0 to 1
  # from http://130.113.54.154/~monger/hsl-rgb.html  
  def _RGB2HSL(self):
    import numpy
    if self.model != 'RGB': return

    HSL = numpy.zeros(3,'d')
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
      HSL[0] = HSL[0]*60.0                                                        # is it necessary to scale to 360 hue values? might be dangerous for small values <1..!
      if (HSL[0] < 0.0):
        HSL[0] = HSL[0] + 360.0
    for i in xrange(2):
      HSL[i+1] = min(HSL[i+1],1.0) 
      HSL[i+1] = max(HSL[i+1],0.0) 
    
    converted = Color('HSL', HSL) 
    self.model = converted.model
    self.color = converted.color
  
  
  # ------------------------------------------------------------------
  # convert R(ed) G(reen) B(lue) to CIE XYZ
  # with all values in the range of 0 to 1
  # from http://www.cs.rit.edu/~ncs/color/t_convert.html
  def _RGB2XYZ(self):
    import numpy
    if self.model != 'RGB': return

    XYZ     = numpy.zeros(3,'d') 
    RGB_lin = numpy.zeros(3,'d')
    convert = numpy.array([[0.412453,0.357580,0.180423],
                           [0.212671,0.715160,0.072169],
                           [0.019334,0.119193,0.950227]])
 
    for i in xrange(3):
      if (self.color[i] > 0.04045): RGB_lin[i] = ((self.color[i]+0.0555)/1.0555)**2.4
      else:                         RGB_lin[i] =   self.color[i]        /12.92
    XYZ = numpy.dot(convert,RGB_lin)
    for i in xrange(3):
      
      XYZ[i] = max(XYZ[i],0.0) 

    converted = Color('XYZ', XYZ) 
    self.model = converted.model
    self.color = converted.color

  
  # ------------------------------------------------------------------
  # convert  CIE XYZ to R(ed) G(reen) B(lue)
  # with all values in the range of 0 to 1
  # from http://www.cs.rit.edu/~ncs/color/t_convert.html
  def _XYZ2RGB(self):
    import numpy
    if self.model != 'XYZ': return

    convert = numpy.array([[ 3.240479,-1.537150,-0.498535],
                           [-0.969256, 1.875992, 0.041556],
                           [ 0.055648,-0.204043, 1.057311]])
    RGB_lin = numpy.dot(convert,self.color)
    RGB     = numpy.zeros(3,'d')

    for i in xrange(3):
      if (RGB_lin[i] > 0.0031308): RGB[i] = ((RGB_lin[i])**(1.0/2.4))*1.0555-0.0555
      else:                        RGB[i] =   RGB_lin[i]             *12.92
    for i in xrange(3):
      RGB[i] = min(RGB[i],1.0) 
      RGB[i] = max(RGB[i],0.0) 
        
    maxVal = max(RGB)                                                                                                                                    # clipping colors according to the display gamut
    if (maxVal > 1.0): RGB /= maxVal
        
    converted = Color('RGB', RGB) 
    self.model = converted.model
    self.color = converted.color

  
  # ------------------------------------------------------------------
  # convert  CIE Lab to CIE XYZ
  # with XYZ in the range of 0 to 1
  # from http://www.easyrgb.com/index.php?X=MATH&H=07#text7
  def _CIELAB2XYZ(self):
    import numpy
    if self.model != 'CIELAB': return
    
    ref_white = numpy.array([.95047, 1.00000, 1.08883],'d')                                                                                              # Observer = 2, Illuminant = D65
    XYZ       = numpy.zeros(3,'d') 

    XYZ[1] = (self.color[0] + 16.0 ) / 116.0
    XYZ[0] = XYZ[1]   + self.color[1]/ 500.0
    XYZ[2] = XYZ[1]   - self.color[2]/ 200.0
    
    for i in xrange(len(XYZ)):
      if (XYZ[i] > 6./29. ): XYZ[i] = XYZ[i]**3.
      else:                  XYZ[i] = 108./841. * (XYZ[i] - 4./29.)
        
    converted = Color('XYZ', XYZ*ref_white) 
    self.model = converted.model
    self.color = converted.color

  
  # ------------------------------------------------------------------    
  # convert CIE XYZ to CIE Lab 
  # with XYZ in the range of 0 to 1
  # from http://en.wikipedia.org/wiki/Lab_color_space, http://www.cs.rit.edu/~ncs/color/t_convert.html
  def _XYZ2CIELAB(self):
    import numpy
    if self.model != 'XYZ': return
    
    ref_white = numpy.array([.95047, 1.00000, 1.08883],'d')                                                                                              # Observer = 2, Illuminant = D65
    XYZ = self.color/ref_white
      
    for i in xrange(len(XYZ)):
      if (XYZ[i] > 216./24389 ): XYZ[i] = XYZ[i]**(1.0/3.0)
      else:                      XYZ[i] = (841./108. * XYZ[i]) + 16.0/116.0
        
    converted = Color('CIELAB', numpy.array([ 116.0 *  XYZ[1] - 16.0,
                                              500.0 * (XYZ[0] - XYZ[1]),
                                              200.0 * (XYZ[1] - XYZ[2]) ]))
    self.model = converted.model
    self.color = converted.color

  
  # ------------------------------------------------------------------                                     
  # convert CIE Lab to Msh colorspace  
  # from http://www.cs.unm.edu/~kmorel/documents/ColorMaps/DivergingColorMapWorkshop.xls
  def _CIELAB2MSH(self):
    import numpy, math
    if self.model != 'CIELAB': return
    
    Msh = numpy.zeros(3,'d') 
    Msh[0] = math.sqrt(numpy.dot(self.color,self.color))
    if (Msh[0] > 0.001):
      Msh[1] = math.acos(self.color[0]/Msh[0])
      if (self.color[1] != 0.0):
        Msh[2] = math.atan2(self.color[2],self.color[1])

    converted = Color('MSH', Msh)
    self.model = converted.model
    self.color = converted.color

  
  # ------------------------------------------------------------------
  # convert Msh colorspace to CIE Lab 
  # s,h in radians
  # from http://www.cs.unm.edu/~kmorel/documents/ColorMaps/DivergingColorMapWorkshop.xls
  def _MSH2CIELAB(self):
    import numpy, math
    if self.model != 'MSH': return
    
    Lab = numpy.zeros(3,'d') 
    Lab[0] = self.color[0] * math.cos(self.color[1])
    Lab[1] = self.color[0] * math.sin(self.color[1]) * math.cos(self.color[2])
    Lab[2] = self.color[0] * math.sin(self.color[1]) * math.sin(self.color[2])

    converted = Color('CIELAB', Lab)
    self.model = converted.model
    self.color = converted.color



### --- COLORMAP CLASS -----------------------------------------------

class Colormap():
  '''
  perceptually uniform diverging or sequential colormaps.
  ''' 

  __slots__ = [
               'left',
               'right',
               'interpolate',
              ]  
  __predefined__ = {
                     'gray':          {'left': Color('HSL',[0,1,1]),          'right': Color('HSL',[0,0,0.15]),          'interpolate': 'perceptualuniform'},
                     'grey':          {'left': Color('HSL',[0,1,1]),          'right': Color('HSL',[0,0,0.15]),          'interpolate': 'perceptualuniform'},
                     'red':           {'left': Color('HSL',[0,1,0.14]),       'right': Color('HSL',[0,0.35,0.91]),       'interpolate': 'perceptualuniform'},
                     'green':         {'left': Color('HSL',[0.33333,1,0.14]), 'right': Color('HSL',[0.33333,0.35,0.91]), 'interpolate': 'perceptualuniform'},
                     'blue':          {'left': Color('HSL',[0.66,1,0.14]),    'right': Color('HSL',[0.66,0.35,0.91]),    'interpolate': 'perceptualuniform'},
                     'seaweed':       {'left': Color('HSL',[0.78,1.0,0.1]),   'right': Color('HSL',[0.40000,0.1,0.9]),   'interpolate': 'perceptualuniform'},
                     'bluebrown':     {'left': Color('HSL',[0.65,0.53,0.49]), 'right': Color('HSL',[0.11,0.75,0.38]),    'interpolate': 'perceptualuniform'},
                     'redgreen':      {'left': Color('HSL',[0.97,0.96,0.36]), 'right': Color('HSL',[0.33333,1.0,0.14]),  'interpolate': 'perceptualuniform'},
                     'bluered':       {'left': Color('HSL',[0.65,0.53,0.49]), 'right': Color('HSL',[0.97,0.96,0.36]),    'interpolate': 'perceptualuniform'},
                     'blueredrainbow':{'left': Color('HSL',[2.0/3.0,1,0.5]),  'right': Color('HSL',[0,1,0.5]),           'interpolate': 'linear'           },
                    }
 
  
  # ------------------------------------------------------------------
  def __init__(self,
               left  = Color('RGB',[1,1,1]),
               right = Color('RGB',[0,0,0]),
               interpolate = 'perceptualuniform',
               predefined = None
               ):

    if str(predefined).lower() in self.__predefined__:
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
    return 'Left: %s Right: %s'%(self.left,self.right)
  
  
  # ------------------------------------------------------------------
  def invert(self):
    (self.left, self.right) = (self.right, self.left)
    return self
  
  
  # ------------------------------------------------------------------  
  def export(self,name = 'uniformPerceptualColorMap',\
                  format = 'paraview',\
                  steps = 2,\
                  crop = [-1.0,1.0],
                  model = 'RGB'):
    ''' 
    [RGB] colormap for use in paraview or gmsh, or as raw string, or array.
    arguments: name, format, steps, crop.
    format is one of (paraview, gmsh, raw, list).
    crop selects a (sub)range in [-1.0,1.0].
    generates sequential map if one limiting color is either white or black,
    diverging map otherwise.
    '''
    
    import copy,numpy,math
    
    def interpolate_Msh(lo, hi, frac):

      def rad_diff(a,b):
        return abs(a[2]-b[2])

      def adjust_hue(Msh_sat, Msh_unsat):                                                                                              # if saturation of one of the two colors is too less than the other, hue of the less
        if Msh_sat[0] >= Msh_unsat[0]:
          return Msh_sat[2]
        else:
          hSpin = Msh_sat[1]/math.sin(Msh_sat[1])*math.sqrt(Msh_unsat[0]**2.0-Msh_sat[0]**2)/Msh_sat[0]
          if Msh_sat[2] < - math.pi/3.0: hSpin *= -1.0
          return Msh_sat[2] + hSpin
     
      Msh1 = numpy.array(lo[:])
      Msh2 = numpy.array(hi[:])

      if (Msh1[1] > 0.05 and Msh2[1] > 0.05 and rad_diff(Msh1,Msh2) > math.pi/3.0): 
        M_mid = max(Msh1[0],Msh2[0],88.0)
        if frac < 0.5:
          Msh2 = numpy.array([M_mid,0.0,0.0],'d')
          frac *= 2.0
        else:
          Msh1 = numpy.array([M_mid,0.0,0.0],'d')
          frac = 2.0*frac - 1.0
      if   Msh1[1] < 0.05 and Msh2[1] > 0.05: Msh1[2] = adjust_hue(Msh2,Msh1)
      elif Msh1[1] > 0.05 and Msh2[1] < 0.05: Msh2[2] = adjust_hue(Msh1,Msh2)
      Msh = (1.0 - frac) * Msh1 + frac * Msh2

      return Color('MSH',Msh)

    def interpolate_linear(lo, hi, frac):
      '''
      linearly interpolate color at given fraction between lower and higher color in model of lower color
      '''

      interpolation = (1.0 - frac) * numpy.array(lo.color[:]) \
                           + frac  * numpy.array(hi.expressAs(lo.model).color[:])

      return Color(lo.model,interpolation)

    def write_paraview(RGB_vector):
      colormap ='<ColorMap name="'+str(name)+'" space="Diverging">\n'
      for i in range(len(RGB_vector)):
        colormap+='<Point x="'+str(i)+'" o="1" r="'+str(RGB_vector[i][0])+'" g="'+str(RGB_vector[i][1])+'" b="'+str(RGB_vector[i][2])+'"/>\n'
      colormap+='</ColorMap>'
      return colormap
    
    def write_gmsh(RGB_vector):
      return 'View.ColorTable = {\n' \
            + ',\n'.join(['{%s}'%(','.join(map(lambda x:str(x*255.0),v))) for v in RGB_vector]) \
            + '\n}'
    
    def write_raw(RGB_vector):
      return '\n'.join(['%s'%('\t'.join(map(lambda x:str(x),v))) for v in RGB_vector])
    
    def write_GOM(RGB_vector):
      return '1 1 '+str(name)+' 9 '+str(name)+' 0 1 0 3 0 0 -1 9 \ 0 0 0 255 255 255 0 0 255 '\
                  + '30 NO_UNIT 1 1 64 64 64 255 1 0 0 0 0 0 0 3 0 ' + str(len(RGB_vector))+' '\
                  .join([' 0 %s 255 1'%(' '.join(map(lambda x:str(int(x*255.0)),v))) for v in reversed(RGB_vector)])+'  '
      
      
    colors = []
    frac = (numpy.array(crop) + 1.0)/2.0

    if self.interpolate == 'perceptualuniform':
      for i in range(steps):
        colors.append(interpolate_Msh(self.left.expressAs('MSH').color,
                                      self.right.expressAs('MSH').color,
                                      float(i)/(steps-1)*(frac[1]-frac[0])+frac[0]))
    elif self.interpolate == 'linear':
      for i in range(steps):
        colors.append(interpolate_linear(self.left,
                                         self.right,
                                         float(i)/(steps-1)*(frac[1]-frac[0])+frac[0]))
    else:
      raise NameError('unknown interpolation method')

    return {\
            'paraview': write_paraview,
            'gmsh':     write_gmsh,
            'gom':      write_GOM,
            'raw':      write_raw,
            'list':     lambda x: x,
    }[format.lower()](map(lambda x:x.expressAs(model).color,colors))
