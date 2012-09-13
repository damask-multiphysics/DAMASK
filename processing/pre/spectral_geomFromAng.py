#!/usr/bin/env python

import os,sys,math,string,numpy
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
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Converts EBSD data from cubic ang files into description for spectral solver (*.geom + material.config)
Can discriminate two phases depending on threshold value 

""" + string.replace('$Id$','\n','\\n')
)


parser.add_option('-c','--column',     dest='column', type='int', \
                                       help='data column to separate phase 1 and 2 [%default]')
parser.add_option('-t','--threshold',  dest='threshold', type='float', \
                                       help='threshold used to separate phases. value < threshold: phase = 1')

parser.set_defaults(column = 1)
parser.set_defaults(threshold = sys.maxint)

(options,filenames) = parser.parse_args()

# ------------------------------------------ setup file handles ---------------------------------------  
eulers = numpy.array([0.0,0.0,0.0],'f')
geomdim = numpy.array([0.0,0.0,0.0],'f')
res = numpy.array([0,0,1],'i')      

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'material':sys.stdout, 'geom':sys.stdout})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append(   {'name':name, 'input':open(name),\
                                    'material':open(os.path.splitext(name)[0]+'_material.config', 'w'),\
                                    'geom':open(os.path.splitext(name)[0]+'.geom', 'w')})

# ------------------------------------------ loop over input files ---------------------------------------  
for file in files:
  point=0
  if file['name'] != 'STDIN': print file['name']
  file['material'].write('#---\n<homogenization>\n#---\n'+
                         '\n[SX]\ntype isostrain\nNgrains 1\n\n'+
                         '#---\n<microstructure>\n#---\n\n')
  tempPart2=             '#---\n<texture>\n---\n\n'
  for line in file['input']:
    lineSplit=line.split()

    if line[0]=='#':
      if len(lineSplit)>2:
        if line.split()[2]=='HexGrid': 
          print 'The file is a hex grid file. Convert it first to sqr grid'
          sys.exit()
        if lineSplit[1]=='XSTEP:':      stepSizeX = float(lineSplit[2]) 
        if lineSplit[1]=='YSTEP:':      stepSizeY = float(lineSplit[2]) 
        if lineSplit[1]=='NCOLS_ODD:':  res[0]    = int(lineSplit[2]) 
        if lineSplit[1]=='NROWS:':      res[1]    = int(lineSplit[2]) 
    else:
      point+=1
      eulers = (float(lineSplit[0])/2.0/math.pi*360.0, \
                   float(lineSplit[1])/2.0/math.pi*360.0, \
                   float(lineSplit[2])/2.0/math.pi*360.0)
               
      if float(lineSplit[options.column-1])<options.threshold:
        phase=1
      else:
        phase=2
        
      file['material'].write(\
        '[Grain%08d]\ncrystallite 1\n(constituent) phase %1d texture %08d fraction 1.0\n' \
                  %(point,phase,point))
      tempPart2+=\
        '[Grain%08d]\n(gauss) phi1 %4.2f Phi %4.2f phi2 %4.2f scatter 0.0 fraction 1.0\n'\
                  %(point,eulers[0],eulers[1],eulers[2])
  geomdim[0] = stepSizeX*res[0]
  geomdim[1] = stepSizeY*res[1]
  geomdim[2] = min(stepSizeX,stepSizeY)
    
  file['material'].write(tempPart2)
  file['geom'].write(
    '3 header\nresolution a %4d b %4d c %1d \ndimension x %5.3f y %5.3f z %5.3f\nhomogenization 1\n'\
               %(res[0],res[1],res[2],geomdim[0],geomdim[1],geomdim[2]))
  file['geom'].write('1 to %d\n'%(res[0]*res[1]))
