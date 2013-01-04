#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

import math

# from http://code.activestate.com/recipes/121574-matrix-vector-multiplication/
def matmult(m, v):
    nrows = len(m)
    w = [None] * nrows
    for row in range(nrows):
        w[row] = reduce(lambda x,y: x+y, map(lambda x,y: x*y, m[row], v))
    return w

# convert H(ue) S(aturation) L(uminance) to R(ot) G(elb) B(lau) 
# with S,L,R,G,B running from 0 to 1, H running from 0 to 360
# from http://en.wikipedia.org/wiki/HSL_and_HSV
def HSL2RGB(HSL):
    RGB = [0.0,0.0,0.0]
    H_strich = HSL[0]/60.0
    c = (1.0- abs(2.0 * HSL[2] - 1.0))*HSL[1]
    x = c*(1.0- abs(H_strich%2-1.0))
    m = HSL[2] -.5*c
    if (0.0 <= H_strich)and(H_strich<1.0):
        RGB[0] = c + m
        RGB[1] = x + m
        RGB[2] = 0.0 + m
    elif (1.0 <= H_strich)and(H_strich<2.0):
        RGB[0] = x + m
        RGB[1] = c + m
        RGB[2] = 0.0 + m
    elif (2.0 <= H_strich)and(H_strich<3.0):
        RGB[0] = 0.0 + m
        RGB[1] = c + m
        RGB[2] = x + m
    elif (3.0 <= H_strich)and(H_strich<4.0):
        RGB[0] = 0.0 + m
        RGB[1] = x + m
        RGB[2] = c + m
    elif (4.0 <= H_strich)and(H_strich<5.0):
        RGB[0] = x + m
        RGB[1] = 0.0 + m
        RGB[2] = c + m
    elif (5.0 <= H_strich)and(H_strich<=6.0):
        RGB[0] = c + m
        RGB[1] = 0.0 + m
        RGB[2] = x + m
    for i in range(3):
        RGB[i] = min(RGB[i],1.0) 
        RGB[i] = max(RGB[i],0.0) 
    return RGB

# convert R(ot) G(elb) B(lau) to H(ue) S(aturation) L(uminance)
# with S,L,R,G,B running from 0 to 1, H running from 0 to 360
# from http://130.113.54.154/~monger/hsl-rgb.html
def RGB2HSL(RGB):
    HSL = [0.0,0.0,0.0]
    maxcolor = max(RGB)
    mincolor = min(RGB)
    HSL[2] = (maxcolor + mincolor)/2.0
    if(mincolor == maxcolor):
        HSL[0] = 0.0
        HSL[1] = 0.0
    else:
        if (HSL[2]<0.5):
            HSL[1] = (maxcolor - mincolor)/(maxcolor + mincolor)
        else:
            HSL[1] = (maxcolor - mincolor)/(2.0 -maxcolor -mincolor)
        if (maxcolor == RGB[0]):
            HSL[0] = 0.0 + (RGB[1] - RGB[2])/(maxcolor - mincolor)
        elif (maxcolor == RGB[1]):
            HSL[0] = 2.0 + (RGB[2] - RGB[0])/(maxcolor - mincolor)
        elif (maxcolor == RGB[2]):
            HSL[0] = 4.0 + (RGB[0] - RGB[1])/(maxcolor - mincolor)
        HSL[0] = HSL[0]*60.0
        if (HSL[0] < 0.0):
            HSL[0] = HSL[0] + 360.0
    for i in range(2):
        HSL[i+1] = min(HSL[i+1],1.0) 
        HSL[i+1] = max(HSL[i+1],0.0) 
    return HSL

# convert R(ot) G(elb) B(lau) to CIE XYZ
# with all values in the range of 0 to 1
# from http://www.cs.rit.edu/~ncs/color/t_convert.html
def RGB2XYZ(RGB):
    XYZ = [0.0,0.0,0.0]
    RGB_lin = [0.0,0.0,0.0]
    for i in range(3):
        if (RGB[i] > 0.04045):
            RGB_lin[i] = ((RGB[i]+0.0555)/1.0555)**2.4
        else:
            RGB_lin[i] = RGB[i]/12.92
    convert =[[0.412453,0.357580,0.180423],[0.212671,0.715160,0.072169],[0.019334,0.119193,0.950227]]
    XYZ = matmult(convert,RGB_lin)
    for i in range(3):
        XYZ[i] = min(XYZ[i],1.0) 
        XYZ[i] = max(XYZ[i],0.0) 
    return XYZ

# convert  CIE XYZ R(ot) G(elb) B(lau)
# with all values in the range of 0 to 1
# from http://www.cs.rit.edu/~ncs/color/t_convert.html
def XYZ2RGB(XYZ):
    RGB_lin = [0.0,0.0,0.0]
    RGB = [0.0,0.0,0.0]
    convert =[[3.240479,-1.537150,-0.498535],[-0.969256,1.875992,0.041556],[0.055648,-0.204043,1.057311]]
    RGB_lin = matmult(convert,XYZ)
    for i in range(3):
        if (RGB_lin[i] > 0.0031308):
            RGB[i] = ((RGB_lin[i])**(1.0/2.4))*1.0555-0.0555
        else:
            RGB[i] = RGB_lin[i]*12.92
    for i in range(3):
        RGB[i] = min(RGB[i],1.0) 
        RGB[i] = max(RGB[i],0.0) 
        
    maxVal = RGB[0]
    if (maxVal < RGB[1]): maxVal = RGB[1]
    if (maxVal < RGB[2]): maxVal = RGB[2]
    if (maxVal > 1.0):
        RGB[0] = RGB[0]/maxVal
        RGB[1] = RGB[1]/maxVal
        RGB[2] = RGB[2]/maxVal
        
    return RGB
    
# convert  CIE Lab to CIE XYZ
# with XYZ in the range of 0 to 1
# from http://www.easyrgb.com/index.php?X=MATH&H=07#text7
def CIELab2XYZ(Lab,white):
    # ref_white_XYZ = [.95047, 1.00000, 1.08883]        # Observer= 2°, Illuminant= D65
    XYZ = [0.0,0.0,0.0]
    var_Y = ( Lab[0] + 16 ) / 116
    var_X = Lab[1] / 500 + var_Y
    var_Z = var_Y - Lab[2] / 200

    if ( var_Y**3 > 0.008856 ):  
        var_Y = var_Y**3
    else: 
        var_Y = ( var_Y-16/116) / 7.787
        
    if ( var_X**3 > 0.008856 ):  
        var_X = var_X**3
    else: 
        var_X = ( var_X-16/116) / 7.787
        
    if ( var_Z**3 > 0.008856 ):  
        var_Z = var_Z**3
    else: 
        var_Z = ( var_Z-16/116) / 7.787

    XYZ[0] = white[0]*var_X          
    XYZ[1] = white[1]*var_Y     
    XYZ[2] = white[2]*var_Z
    
    return XYZ  
  
# convert CIE XYZ to CIE Lab 
# with XYZ in the range of 0 to 1
# from http://en.wikipedia.org/wiki/Lab_color_space, http://www.cs.rit.edu/~ncs/color/t_convert.html
def XYZ2CIELab(XYZ,white):
    # ref_white_XYZ = [.95047, 1.00000, 1.08883]          # Observer= 2°, Illuminant= D65
    Lab = [0.0,0.0,0.0]
    var_X = XYZ[0]/white[0] 
    var_Y = XYZ[1]/white[1]          
    var_Z = XYZ[2]/white[2] 
    
    if ( var_X > 0.008856 ):
        var_X = var_X**(1.0/3.0)
    else: 
        var_X = ( 7.787 * var_X ) + (16.0/116.0)
        
    if ( var_Y > 0.008856 ):
        var_Y = var_Y**(1.0/3.0)
    else:
        var_Y = ( 7.787 * var_Y ) + (16.0/116.0)
        
    if ( var_Z > 0.008856 ):
        var_Z = var_Z**(1.0/3.0)
    else: 
        var_Z = (7.787*var_Z)+(16.0/116.0)
        
    Lab[0] = (116.0 * var_Y) - 16.0
    Lab[1] = 500.0 * (var_X - var_Y)
    Lab[2] = 200.0 * (var_Y - var_Z)
    return Lab
 
# convert Cie Lab to msh colorspace  
# from http://www.cs.unm.edu/~kmorel/documents/ColorMaps/DivergingColorMapWorkshop.xls
def CIELab2Msh(Lab):
    Msh = [0.0,0.0,0.0]
    Msh[0] = math.sqrt(Lab[0]**2.0 + Lab[1]**2.0 + Lab[2]**2.0)
    if (Msh[0] != 0.0) and (Msh[0] > 0.001):
        Msh[1] = math.acos(Lab[0]/Msh[0])
    if (Lab[1] != 0.0) and (Msh[1] > 0.001):
        Msh[2] = math.atan2(Lab[2],Lab[1])
    return Msh

# convert  msh colorspace to Cie Lab 
# from http://www.cs.unm.edu/~kmorel/documents/ColorMaps/DivergingColorMapWorkshop.xls
def Msh2CIELab(Msh):
    Lab = [0.0,0.0,0.0]
    Lab[0] = Msh[0] * math.cos(Msh[1])
    Lab[1] = Msh[0] * math.sin(Msh[1]) * math.cos(Msh[2])
    Lab[2] = Msh[0] * math.sin(Msh[1]) * math.sin(Msh[2])
    return Lab
