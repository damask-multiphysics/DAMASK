#!/usr/bin/env python

import math, string, sys
from damask import Colormaps
from optparse import OptionParser, Option

# -----------------------------
class extendableOption(Option):
# -----------------------------
# used for definition of new option parser action 'extend', which enables to take multiple option arguments
# taken from online tutorial http://docs.python.org/library/optparse.html
  
  ACTIONS = Option.ACTIONS + ("extend",)
  STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
  TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
  ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

  def take_action(self, action, dest, opt, value, values, parser):
    if action == "extend":
      lvalue = value.split(",")
      values.ensure_value(dest, []).extend(lvalue)
    else:
      Option.take_action(self, action, dest, opt, value, values, parser)



# --------------------------------------------------------------------
                               # MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing Cauchy stress based on given column(s) of
deformation gradient and first Piola--Kirchhoff stress.

""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-f','--farbe', dest='farbe', type='float', nargs=3, \
                  help='RGB values of the desired color')
parser.add_option('-c','--colormodel', dest='colormodel', \
                  help='colormodel of left and right "RGB","RGB255","HSV","HSL" [%default]')
parser.add_option('-o','--outtype', dest='outtype', \
                  help='output file type "paraview","gmsh","raw" [%default]')
parser.add_option('-s','--steps', dest='steps', type='int', nargs = 1, \
                  help='no of interpolation steps [%default]')
parser.add_option('-m','--maptype', dest='maptype', \
                  help='Increasing or decreasing sequential map "inc","dec" [%default]')									

parser.set_defaults(colormodel = 'RGB')
parser.set_defaults(outtype = 'paraview')
parser.set_defaults(steps = '10')
parser.set_defaults(maptype = 'dec')

(options,filenames) = parser.parse_args()


# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
  
  
# -----------------------------------------------------------------------------------------------------  

myColorMap = Colormaps()
   
def rad_dif(Msh1,Msh2,white):
    HSL1 = myColorMap.RGB2HSL(myColorMap.XYZ2RGB(myColorMap.CIELab2XYZ(myColorMap.Msh2CIELab(Msh1),white)))
    HSL2 = myColorMap.RGB2HSL(myColorMap.XYZ2RGB(myColorMap.CIELab2XYZ(myColorMap.Msh2CIELab(Msh2),white)))
    return abs(HSL1[0]*math.pi/180.0-HSL2[0]*math.pi/180.0)

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
		Msh_mid = [0.0,0.0,0.0]
		if (options.maptype.lower() == 'dec'):
				Msh1 = myColorMap.CIELab2Msh(myColorMap.XYZ2CIELab(myColorMap.RGB2XYZ(RGB1),white))
				Msh2 = myColorMap.CIELab2Msh(myColorMap.XYZ2CIELab(myColorMap.RGB2XYZ(white),white))
				if (Msh1[1] > 0.05):
						Msh_mid[0] = max(Msh1[0],Msh2[0],88.0)
				if (Msh1[1] < 0.05):
						Msh1[2] = adjust_hue(Msh2,Msh1[0])
		if (options.maptype.lower() == 'inc'):
				Msh1 = myColorMap.CIELab2Msh(myColorMap.XYZ2CIELab(myColorMap.RGB2XYZ(white),white))
				Msh2 = myColorMap.CIELab2Msh(myColorMap.XYZ2CIELab(myColorMap.RGB2XYZ(RGB1),white))
				if (Msh2[1] > 0.05):
						Msh_mid[0] = max(Msh1[0],Msh2[0],88.0)
				if (Msh2[1] < 0.05):
						Msh2[2] = adjust_hue(Msh2,Msh1[0])
		for i in range(3):
				Msh_mid[i] = (1.0-interp)*Msh1[i] + interp* Msh2[i]
		return myColorMap.XYZ2RGB(myColorMap.CIELab2XYZ(myColorMap.Msh2CIELab(Msh_mid),white))

white = [0.950456, 1.0, 1.088754]
interpolatorArray = []
for i in range(options.steps+1): interpolatorArray.append(float(i)/options.steps)
rMatrix = []
gMatrix = []
bMatrix = []
for i in interpolatorArray:
    step_no = str(interpolatorArray.index(i))
    color = interpolate_color(options.farbe,white,i)
    rMatrix.append(color[0])
    gMatrix.append(color[1])
    bMatrix.append(color[2])

colorMatrix = [rMatrix,gMatrix,bMatrix]    

if options.outtype.lower() == 'paraview':
    myColorMap.write_paraview(colorMatrix,filenames[0])

