#!/usr/bin/env python

import math, string, sys
from damask import Color,Colormap
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

parser = OptionParser(option_class=extendableOption, usage='%prog options [file]', description = """
Generate uniform perceptual colormap.

""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-l','--left', dest='left', type='float', nargs=3, \
                  help='left color [%default]')
parser.add_option('-r','--right', dest='right', type='float', nargs=3, \
                  help='right color [%default]')
parser.add_option('-c','--colormodel', dest='colormodel', \
                  help='colormodel of left and right "RGB","HSL","XYZ","CIELAB","MSH" [%default]')
parser.add_option('-f','--format', dest='format', \
                  help='output file format "paraview","gmsh","raw" [%default]')
parser.add_option('-s','--steps', dest='steps', type='int', nargs = 1, \
                  help='no of interpolation steps [%default]')
parser.add_option('-t','--trim', dest='trim', type='float', nargs = 2, \
                  help='trim the colormap w.r.t the given values [%default]')

parser.set_defaults(colormodel = 'RGB')
parser.set_defaults(format = 'paraview')
parser.set_defaults(steps = 10)
parser.set_defaults(trim = [-1.0,1.0])
parser.set_defaults(left = [1.0,1.0,1.0])
parser.set_defaults(right = [0.0,0.0,0.0])
(options,filenames) = parser.parse_args()

  
  
# -----------------------------------------------------------------------------------------------------  

leftColor = Color(options.colormodel.upper(),list(options.left))
rightColor = Color(options.colormodel.upper(),list(options.right))

myColormap = Colormap(leftColor,rightColor)


if filenames == []:
  outFile = sys.stdout
  outColormap = myColormap.export('uniform perceptual colormap',options.format.lower(),options.steps,list(options.trim))
else:
  outColormap = myColormap.export(filenames[0],options.format.lower(),options.steps,list(options.trim))
  if options.format.lower() == 'paraview':
    if filenames[0].endswith('.xml'):
      outFile = open('%s'%filenames[0],'w')
    else:
      outFile = open('%s.xml'%filenames[0],'w')

  elif options.format.lower() == 'gmsh':
    if filenames[0].endswith('.msh'): 
      outFile = open('%s'%filenames[0],'w')
    else:
      outFile = open('%s.msh'%filenames[0],'w')

  elif options.format.lower() == 'raw':
    if filenames[0].endswith('.txt'): 
      outFile = open('%s'%filenames[0],'w')
    else:
      outFile = open('%s.txt'%filenames[0],'w')

outFile.write(outColormap)
outFile.close()