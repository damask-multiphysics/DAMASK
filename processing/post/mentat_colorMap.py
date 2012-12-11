#!/usr/bin/env python 

import sys, os, string, damask
from colorsys import *
from optparse import OptionParser

sys.path.append(damask.solver.Marc().libraryPath('../../'))

try:
  from py_mentat import *
except:
  print('warning: no valid Mentat release found')


# -----------------------------
def readConfig(configFile,ownPath):
  config = {}
  configDir = os.path.split(os.path.realpath(ownPath))[0]
  filename = os.path.join(configDir,configFile)
  if os.path.isfile(filename):
    file = open(filename)
    content = file.readlines()
    file.close()
    for line in content:
      item = line.split()
      config[item[0]] = {}
      config[item[0]]['lower'] = map(float,item[1].split(','))
      config[item[0]]['upper'] = map(float,item[2].split(','))
      config[item[0]]['symmetric'] = len(item) > 3
  
  return config



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
def interpolate(val0, val1, x):
    return val0 + (val1 - val0) * x



# -----------------------------
def syminterpolate(comp, val0, val1, x):
    if comp == "hue":
        return {True:val0,False:val1}[x<0.5]

    if comp == "lightness":
        val_middle = 1
    elif comp == "saturation":
        val_middle = 0
    if x < 0.5:
        return interpolate(val0, val_middle, 2*x)
    else:
        return interpolate(val_middle, val1, 2*x-1)
        
        

# -----------------------------
def colorMap(colors,baseIdx=32):
    cmds = [ "*color %i %f %f %f"%(idx+baseIdx,color[0],color[1],color[2]) 
             for idx,color in enumerate(colors) ]    
    return cmds
    

# -----------------------------
# MAIN FUNCTION STARTS HERE
# -----------------------------

parser = OptionParser(usage="%prog [options] configured scheme | (lower_h,l,s upper_h,l,s)", description = """
Changes the color map in mentat. 

Interpolates colors between "lower_hls" and "upper_hls". 
For symmetric scales use option "-s". 

Example colors:
- Non-symmetric scales: 0.167,0.9,0.1  0.167,0.1,0.9
- Symmetric scales: 0,0.2,0.9  0.333,0.2,0.9
""" + string.replace('$Id$','\n','\\n')
)


parser.add_option("-s","--symmetric", action = "store_true", 
                  dest = "symmetric", \
                  help = "symmetric legend [%default]")
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
parser.add_option("-c", "--config", type='string', \
                  dest = "config",\
                  help = "configuration file [%default]")
parser.add_option("-v", "--verbose", action="store_true",\
                  dest = "verbose",\
                  help = "write Mentat command stream also to stdout [%default]")

parser.set_defaults(port = 40007)
parser.set_defaults(baseIdx = 32)
parser.set_defaults(colorcount = 32)
parser.set_defaults(config = 'colorMap.config')
parser.set_defaults(symmetric = False)
parser.set_defaults(inverse   = False)
parser.set_defaults(palette   = False)
parser.set_defaults(verbose   = False)

msg = []

(options, colors) = parser.parse_args()

config = readConfig(options.config,sys.argv[0])

if len(colors) == 0:
  hlsColor_range = (options.symmetric and [[0,0.2,0.9],[0.333,0.2,0.9]]) or \
                                          [[0.167,0.9,0.1],[0.167,0.1,0.9]]
elif len(colors) == 2:
	hlsColor_range = [map(float, colors[i].split(',')) for i in range(2)]
elif colors[0] in config:
  options.symmetric = config[colors[0]]['symmetric']
  hlsColor_range = [config[colors[0]]['lower'],\
  									config[colors[0]]['upper']]
else:
  msg.append('two color tuples required')

if msg == []:  
	hlsColors_limits = [[0.0,0.0,0.0],[1.0,1.0,1.0]]

	if options.inverse:
	  hlsColor_range = [hlsColor_range[1],hlsColor_range[0]]

	for i in range(2):			
		for j in range(min(3,len(hlsColor_range[i]))):
			if hlsColor_range[i][j] < hlsColors_limits[0][j] or hlsColor_range[i][j] > hlsColors_limits[1][j]:
				msg.append('%s of %s color exceeds limit'%(['hue','lightness','saturation'][j],limit))
		
if msg != []:
    parser.error('\n'+'\n'.join(msg)+'\n')

### interpolate hls values

if options.symmetric:
    hlsColors = [ [ syminterpolate(comp, hlsColor_range[0][j], hlsColor_range[1][j], float(idx)/(options.colorcount-1)) 
                    for j,comp in enumerate(["hue","lightness","saturation"]) ]
                  for idx in range(options.colorcount) ]
else:
    hlsColors = [ [ interpolate(hlsColor_range[0][j], hlsColor_range[1][j], float(idx)/(options.colorcount-1)) 
                    for j,comp in enumerate(["hue","lightness","saturation"]) ]
                  for idx in range(options.colorcount) ]



### convert to rgb values

rgbColors = [ hls_to_rgb(hlsColor[0], hlsColor[1], hlsColor[2]) 
              for hlsColor in hlsColors ]

if options.palette:
  for rgb in rgbColors:
    print '\t'.join(map(lambda x: str(int(255*x)),rgb))
  sys.exit(0)
if options.palettef:
  for rgb in rgbColors:
    print '\t'.join(map(str,rgb))
  sys.exit(0)
  
### connect to mentat and change colorMap

outputLocals = {}
print 'waiting to connect...'
py_connect('',options.port)
print 'connected...'

cmds = colorMap(rgbColors,options.baseIdx)
output(['*show_table']+cmds+['*show_model *redraw'],outputLocals,'Mentat')
py_disconnect()

if options.verbose:
  output(cmds,outputLocals,'Stdout')
  
  
