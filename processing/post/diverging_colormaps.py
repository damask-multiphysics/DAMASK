#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

import math, convert_colormodels
def rad_dif(Msh1,Msh2,white):
    HSL1 = convert_colormodels.RGB2HSL(convert_colormodels.XYZ2RGB(convert_colormodels.CIELab2XYZ(convert_colormodels.Msh2CIELab(Msh1),white)))
    HSL2 = convert_colormodels.RGB2HSL(convert_colormodels.XYZ2RGB(convert_colormodels.CIELab2XYZ(convert_colormodels.Msh2CIELab(Msh2),white)))
    return abs(HSL1[0]*math.pi/180.0-HSL2[0]*math.pi/180.0)

def adjust_hue(Msh_sat,M_unsat):
    if Msh_sat[0] >= M_unsat:
        return Msh_sat[2]
    else:
        hSpin = Msh_sat[1]*math.sqrt((M_unsat)**2.0-(Msh_sat[0])**2)/(Msh_sat[0]*math.sin(Msh_sat[1]))
        if Msh_sat[2] > - math.pi/3.0:
            return Msh_sat[2] + hSpin
        else:
            return Msh_sat[2] - hSpin
        
    
def interpolate_color(RGB1,RGB2,white,interp):
    Msh1 = convert_colormodels.CIELab2Msh(convert_colormodels.XYZ2CIELab(convert_colormodels.RGB2XYZ(RGB1),white))
    Msh2 = convert_colormodels.CIELab2Msh(convert_colormodels.XYZ2CIELab(convert_colormodels.RGB2XYZ(RGB2),white))
    Msh_mid = [0.0,0.0,0.0]
    if ((Msh1[1] > 0.05 and Msh2[1] > 0.05) and rad_dif(Msh1,Msh2,white) > math.pi/3.0):
        Msh_mid[0] = max(Msh1[0],Msh2[0],88.0)
        if interp < 0.5:
            Msh2[0] = Msh_mid[0]
            Msh2[1] = 0.0
            Msh2[2] = 0.0
            interp = 2.0*interp
        else:
            Msh1[0] = Msh_mid[0]
            Msh1[1] = 0.0
            Msh1[2] = 0.0
            interp = 2.0*interp - 1.0
    if (Msh1[1] < 0.05) and (Msh2[1] > 0.05):
        Msh1[2] = adjust_hue(Msh2,Msh1[0])
    elif (Msh2[1] < 0.05) and (Msh1[1] > 0.05):
        Msh2[2] = adjust_hue(Msh1,Msh2[0])
    for i in range(3):
        Msh_mid[i] = (1.0-interp)*Msh1[i] + interp* Msh2[i]
    return convert_colormodels.XYZ2RGB(convert_colormodels.CIELab2XYZ(convert_colormodels.Msh2CIELab(Msh_mid),white))

test1 = [0.231372549,0.298039216,0.752941176]
test2 = [0.705882353,0.015686275,0.149019608]
#test1 = [0.0,1.0,24.0/255.0]
#test1 = [0.5,0.5,0.5]
#test2 = [0.0,151.0/255.0,21.0/255.0]
x=0.950456 
y=1.0
z=1.088754
white = [x, y, z]
test = test1
iteration = 33
delta = 1.0/float(iteration-1)
f = -delta
for i in range(iteration):
    f = f + delta
    test = test + interpolate_color(test1,test2,white,f)
test = test + test2
for i in range(iteration):
    print i
    print test[3*i+3]*255.0
    print test[3*i+1+3]*255.0
    print test[3*i+2+3]*255.0, '\n'
