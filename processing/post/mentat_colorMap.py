#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os, sys, string
import damask
from optparse import OptionParser

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# -----------------------------
def outMentat(cmd,locals):
    if cmd[0:3] == '(!)':
        exec(cmd[3:])
    elif cmd[0:3] == '(?)':
        cmd = eval(cmd[3:])
        py_send(cmd)
    else:
        py_send(cmd)
    return



# -----------------------------
def outStdout(cmd,locals):
    if cmd[0:3] == '(!)':
        exec(cmd[3:])
    elif cmd[0:3] == '(?)':
        cmd = eval(cmd[3:])
        print cmd
    else:
        print cmd
    return



# -----------------------------
def output(cmds,locals,dest):
    for cmd in cmds:
        if isinstance(cmd,list):
            output(cmd,locals,dest)
        else:
            {\
            'Mentat': outMentat,\
            'Stdout': outStdout,\
            }[dest](cmd,locals)
    return



# -----------------------------
def colorMap(colors,baseIdx=32):
    cmds = [ "*color %i %f %f %f"%(idx+baseIdx,color[0],color[1],color[2]) 
             for idx,color in enumerate(colors) ]    
    return cmds
    

# -----------------------------
# MAIN FUNCTION STARTS HERE
# -----------------------------

parser = OptionParser(option_class=damask.extendableOption, usage="%prog [options] predefinedScheme | (lower_h,s,l upper_h,s,l)", description = """
Changes the color map in MSC.Mentat. 

Interpolates colors between "lower_hsl" and "upper_hsl". 

""", version = scriptID)

parser.add_option("-i","--inverse", action = "store_true", 
                  dest = "inverse", \
                  help = "invert legend [%default]")
parser.add_option(     "--palette", action = "store_true", 
                  dest = "palette", \
                  help = "output plain rgb palette integer values (0-255) [%default]")
parser.add_option(     "--palettef", action = "store_true", 
                  dest = "palettef", \
                  help = "output plain rgb palette float values (0.0-1.0) [%default]")
parser.add_option("-p", "--port", type = "int",\
                  dest = "port",\
                  help = "Mentat connection port [%default]")
parser.add_option("-b", "--baseindex", type = "int",\
                  dest = "baseIdx",\
                  help = "base index of colormap [%default]")
parser.add_option("-n", "--colorcount", type = "int",\
                  dest = "colorcount",\
                  help = "number of colors [%default]")
parser.add_option("-v", "--verbose", action="store_true",\
                  dest = "verbose",\
                  help = "write Mentat command stream also to STDOUT [%default]")

parser.set_defaults(port = 40007)
parser.set_defaults(baseIdx = 32)
parser.set_defaults(colorcount = 32)
parser.set_defaults(inverse   = False)
parser.set_defaults(palette   = False)
parser.set_defaults(palettef  = False)
parser.set_defaults(verbose   = False)

msg = []

(options, colors) = parser.parse_args()

if len(colors) == 0:
  parser.error('missing color information')
  
elif len(colors) == 1:
  theMap = damask.Colormap(predefined = colors[0])

elif len(colors) == 2:
  theMap = damask.Colormap(damask.Color('HSL',map(float, colors[0].split(','))),
                           damask.Color('HSL',map(float, colors[1].split(','))) )

else:
  theMap = damask.Colormap()

if options.inverse:
  theMap = theMap.invert()

if options.palettef:
  print theMap.export(format='raw',steps=options.colorcount)
elif options.palette:
  for theColor in theMap.export(format='list',steps=options.colorcount):
    print '\t'.join(map(lambda x: str(int(255*x)),theColor))
else:  
### connect to Mentat and change colorMap
  sys.path.append(damask.solver.Marc().libraryPath('../../'))
  try:
    from py_mentat import *
    print 'waiting to connect...'
    py_connect('',options.port)
    print 'connected...'
    mentat = True
  except:
    sys.stderr.write('warning: no valid Mentat release found\n')
    mentat = False

  outputLocals = {}
  cmds = colorMap(theMap.export(format='list',steps=options.colorcount),options.baseIdx)
  if mentat:
    output(['*show_table']+cmds+['*show_model *redraw'],outputLocals,'Mentat')
    py_disconnect()
  
  if options.verbose:
    output(cmds,outputLocals,'Stdout')
