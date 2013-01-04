#!/usr/bin/env python

import math, convert_colormodels, colormap_io, string, sys

    
def adjust_hue(Msh_sat,M_unsat):
    if ( Msh_sat[0] >= (M_unsat-0.1) ):
        return Msh_sat[2]
    else:
        hSpin = Msh_sat[1]*math.sqrt((M_unsat)**2.0-(Msh_sat[0])**2)/(Msh_sat[0]*math.sin(Msh_sat[1]))
        if Msh_sat[2] > - math.pi/3.0:
            return Msh_sat[2] + hSpin
        else:
            return Msh_sat[2] - hSpin


def interpolate_color(RGB1,white,interp):
    Msh1 = convert_colormodels.CIELab2Msh(convert_colormodels.XYZ2CIELab(convert_colormodels.RGB2XYZ(RGB1),white))
    Msh2 = [0.0,0.0,0.0]
    Msh_mid = [0.0,0.0,0.0]
    if ((Msh1[1] > 0.05 and Msh2[1] > 0.05):
        Msh_mid[0] = max(Msh1[0],Msh2[0],88.0)
        if interp <=1:
            Msh2[0] = Msh_mid[0]
            Msh2[1] = 0.0
            Msh2[2] = 0.0
            interp = 2.0*interp
        else:
            print 'Interpolation factor value must be within 0 and 1!'
    if (Msh1[1] < 0.05):
        Msh1[2] = adjust_hue(Msh2,Msh1[0])
    else:
        print 'The Saturation value of the given color is aggreable!'
    for i in range(3):
        Msh_mid[i] = (1.0-interp)*Msh1[i] + interp* Msh2[i]
    return convert_colormodels.XYZ2RGB(convert_colormodels.CIELab2XYZ(convert_colormodels.Msh2CIELab(Msh_mid),white))            
    
    
ex1 = [46/255 139/255 87/255]

interpolatorArray = []
for i in range(options.steps+1): interpolatorArray.append(float(i)/options.steps)
rMatrix = []
gMatrix = []
bMatrix = []
for i in interpolatorArray:
    step_no = str(interpolatorArray.index(i))
    color = interpolate_color(options.left,options.right,white,i)
    rMatrix.append(color[0])
    gMatrix.append(color[1])
    bMatrix.append(color[2])
    print 'step no: %s'%step_no
    print color[0], color[1], color[2]

colorMatrix = [rMatrix,gMatrix,bMatrix]    