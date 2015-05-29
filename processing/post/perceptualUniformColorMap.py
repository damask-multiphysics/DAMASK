#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import math, string, sys, os
import damask
from optparse import OptionParser

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
                               # MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Produces perceptually linear diverging and sequential colormaps in formats suitable for visualization software or simply as a list of interpolated colors.

""", version = scriptID)

parser.add_option('-l','--left', dest='left', type='float', nargs=3, \
                  help='left color [%default]')
parser.add_option('-r','--right', dest='right', type='float', nargs=3, \
                  help='right color [%default]')
parser.add_option('-c','--colormodel', dest='colormodel', \
                  help='colormodel of left and right "RGB","HSL","XYZ","CIELAB","MSH" [%default]')
parser.add_option('-p','--predefined', dest='predefined', \
                  help='predefined colormap [%default]')
parser.add_option('-f','--format', dest='format', action='extend', \
                  help='output file format "paraview","gmsh","raw","GOM",[paraview, autodetect if output file extension is given]')
parser.add_option('-s','--steps', dest='steps', type='int', nargs = 1, \
                  help='number of interpolation steps [%default]')
parser.add_option('-t','--trim', dest='trim', type='float', nargs = 2, \
                  help='trim the colormap w.r.t the given values [%default]')
parser.set_defaults(colormodel = 'RGB')
parser.set_defaults(predefined = None)
parser.set_defaults(format = [''])
parser.set_defaults(steps = 10)
parser.set_defaults(trim  = (-1.0,1.0))
parser.set_defaults(left  = (1.0,1.0,1.0))
parser.set_defaults(right = (0.0,0.0,0.0))
(options,filenames) = parser.parse_args()

outtypes   = ['paraview','gmsh','raw','GOM']
extensions = ['.xml','.msh','.txt','.legend']
if options.trim[0]< -1.0 or \
   options.trim[1] > 1.0 or \
   options.trim[0]>= options.trim[1]:
  parser.error('invalid trim range ...')

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == [] and options.format == ['']:
  files.append({'outtype':'paraview','output':sys.stdout,'name':'colormap'})

if (len(options.format) == (len(filenames)+1)) and (len(options.format) > 1):
  for i in xrange(1,len(options.format)):
    [basename,myExtension] = os.path.splitext(os.path.basename(filenames[i-1]))
    if options.format[i] in outtypes:
      myExtension = extensions[outtypes.index(options.format[i])]
      myType = outtypes[extensions.index(myExtension)]
    elif myExtension in extensions:
      myType = outtypes[extensions.index(myExtension)]
    else:
      myType = 'paraview'
      myExtension = extensions[outtypes.index(myType)]      
    files.append({'name':basename, 'output':open(basename+myExtension,'w'), 'outtype': myType})

if (len(options.format) > (len(filenames)+1)): 
  if (len(filenames) == 1) :
    [basename,myExtension] = os.path.splitext(os.path.basename(filenames[0]))
    for i in xrange(1,len(options.format)):
      if options.format[i] in outtypes:
        myExtension = extensions[outtypes.index(options.format[i])]
        myType = outtypes[extensions.index(myExtension)]
      else:
        myType = 'paraview'
        myExtension = extensions[outtypes.index(myType)]
      files.append({'name':basename, 'output':open(basename+myExtension,'w'), 'outtype': myType})
  elif len(filenames) == 0:
    for i in xrange(1,len(options.format)):
      if options.format[i] in outtypes:
        myExtension = extensions[outtypes.index(options.format[i])]
        myType = outtypes[extensions.index(myExtension)]
        basename = myType
      else:
        myType = 'paraview'
        myExtension = extensions[outtypes.index(myType)]
        basename = myType    
      files.append({'name':basename, 'output':open(basename+myExtension,'w'), 'outtype': myType})
  elif len(filenames) > 1:
    for i in xrange(len(filenames)):
      [basename,myExtension] = os.path.splitext(os.path.basename(filenames[i]))
      if options.format[i+1] in outtypes:
        myExtension = extensions[outtypes.index(options.format[i+1])]
        myType = outtypes[extensions.index(myExtension)]
      elif myExtension in extensions:
        myType = outtypes[extensions.index(myExtension)]
      else:
        myType = 'paraview'
        myExtension = extensions[outtypes.index(myType)]
      files.append({'name':basename, 'output':open(basename+myExtension,'w'), 'outtype': myType})
    for i in xrange(len(filenames)+1,len(options.format)):
      if options.format[i] in outtypes:
        myExtension = extensions[outtypes.index(options.format[i])]
        myType = outtypes[extensions.index(myExtension)]
        basename = myType
      else:
        myType = 'paraview'
        myExtension = extensions[outtypes.index(myType)]
        basename = myType
      files.append({'name':basename, 'output':open(basename+myExtension,'w'), 'outtype': myType})

if (len(options.format) < (len(filenames)+1)) and (options.format!=['']):
  for i in xrange(1,len(options.format)):
    [basename,myExtension] = os.path.splitext(os.path.basename(filenames[i-1]))
    if options.format[i] in outtypes:
      myExtension = extensions[outtypes.index(options.format[i])]
      myType = outtypes[extensions.index(myExtension)]
    elif (options.format[i] not in outtypes) and (myExtension in extensions):
      myType = outtypes[extensions.index(myExtension)]
    else:
      myType = 'paraview'
      myExtension = extensions[outtypes.index(myType)]
    files.append({'name':basename, 'output':open(basename+myExtension,'w'), 'outtype': myType})      # files.append({'name':basename, 'output':open(basename+myExtension,'w'), 'outtype': myType})
  for i in xrange((len(options.format)-1),len(filenames)):
    [basename,myExtension] = os.path.splitext(os.path.basename(filenames[i]))
    if myExtension.lower() in extensions:
      myType = outtypes[extensions.index(myExtension)]
    else:
      myType = 'paraview'
      myExtension = extensions[outtypes.index(myType)]
    files.append({'name':basename, 'output':open(basename+myExtension,'w'), 'outtype': myType})
    
elif (len(filenames)> 0) and (options.format == ['']):
 for [i,name] in enumerate(filenames):
  [basename,myExtension] = os.path.splitext(os.path.basename(name))
  if myExtension.lower() in extensions:
    myType = outtypes[extensions.index(myExtension)]
  else:
    myType = 'paraview'
    myExtension = extensions[outtypes.index(myType)]
  files.append({'name':basename, 'output':open(basename+myExtension,'w'), 'outtype': myType})  
    
leftColor = damask.Color(options.colormodel.upper(),list(options.left))
rightColor = damask.Color(options.colormodel.upper(),list(options.right))
myColormap = damask.Colormap(leftColor,rightColor,predefined = options.predefined)

for file in files:
  outColormap = myColormap.export(file['name'],file['outtype'],options.steps,list(options.trim))
  file['output'].write(outColormap)
  file['output'].close()  
