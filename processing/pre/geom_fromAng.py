#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math,string
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """

Generate geometry description and material configuration from EBSD data in given square-gridded 'ang' file.
Two phases can be discriminated based on threshold value in a given data column.

""", version = scriptID)

parser.add_option('--column',              dest='column', type='int', metavar = 'int', \
                  help='data column to discriminate between both phases [%default]')
parser.add_option('-t','--threshold',      dest='threshold', type='float', metavar = 'float', \
                  help='threshold value for phase discrimination [%default]')
parser.add_option('--homogenization',      dest='homogenization', type='int', metavar = 'int', \
                  help='homogenization index for <microstructure> configuration [%default]')
parser.add_option('--phase',               dest='phase', type='int', nargs = 2, metavar = 'int int', \
                  help='phase indices for <microstructure> configuration %default')
parser.add_option('--crystallite',         dest='crystallite', type='int', metavar = 'int', \
                  help='crystallite index for <microstructure> configuration [%default]')
parser.add_option('-c', '--configuration', dest='config', action='store_true', \
                  help='output material configuration [%default]')
parser.add_option('--compress',            dest='compress', action='store_true', \
                  help='lump identical microstructure and texture information [%default]')
parser.add_option('-a', '--axes',         dest='axes', nargs = 3, metavar = 'string string string', \
                  help='Euler angle coordinate system for <texture> configuration x,y,z = %default')
    
                  
parser.set_defaults(column         = 11)
parser.set_defaults(threshold      = 0.5)
parser.set_defaults(homogenization = 1)
parser.set_defaults(phase          = [1,2])
parser.set_defaults(crystallite    = 1)
parser.set_defaults(config         = False)
parser.set_defaults(compress       = False)
parser.set_defaults(axes           = ['y','x','-z'])
(options,filenames) = parser.parse_args()

for i in options.axes:
  if i.lower() not in ['x','+x','-x','y','+y','-y','z','+z','-z']:
    parser.error('invalid axes %s %s %s' %(options.axes[0],options.axes[1],options.axes[2]))

#--- setup file handles ---------------------------------------------------------------------------
files = []
if filenames == []:
  files.append({'name':'STDIN',
                'input':sys.stdin,
                'output':sys.stdout,
                'croak':sys.stderr,
               })
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name,
                    'input':open(name),
                    'output':open(name+'_tmp','w'),
                    'croak':sys.stdout,
                    })


#--- loop over input files ------------------------------------------------------------------------
for file in files:
  file['croak'].write('\033[1m' + scriptName + '\033[0m: ' + (file['name'] if file['name'] != 'STDIN' else '') + '\n')

  info = {
          'grid':   np.ones (3,'i'),
          'size':   np.zeros(3,'d'),
          'origin': np.zeros(3,'d'),
          'microstructures': 0,
          'homogenization':  options.homogenization,
         }

  step  = [0,0]
  point = 0
  for line in file['input']:
    words = line.split()
    if len(words) == 0: continue                                                                    # ignore empty lines
    if words[0] == '#':                                                                             # process initial comments block
      if len(words) > 2:
        if words[2].lower() == 'hexgrid': 
          file['croak'].write('The file has HexGrid format. Please first convert to SquareGrid...\n')
          break
        if words[1] == 'XSTEP:':     step[0] = float(words[2])
        if words[1] == 'YSTEP:':     step[1] = float(words[2])
        if words[1] == 'NCOLS_ODD:':
          info['grid'][0] = int(words[2])
          eulerangles = np.zeros((info['grid'][0]*info['grid'][1],3),dtype='f')
          phase       = np.zeros(info['grid'][0]*info['grid'][1],dtype='i')
        if words[1] == 'NROWS:':
          info['grid'][1] = int(words[2])
          eulerangles = np.zeros((info['grid'][0]*info['grid'][1],3),dtype='f')
          phase       = np.zeros(info['grid'][0]*info['grid'][1],dtype='i')
    else:                                                                                           # finished with comments block
      phase[point] = options.phase[int(float(words[options.column-1]) > options.threshold)]
      eulerangles[point,...] = map(lambda x: math.degrees(float(x)), words[:3])
      point += 1
 
  if info['grid'].prod() != point:
    file['croak'].write('Error: found %s microstructures. Header info in ang file might be wrong.\n'%point)
    continue
  if options.compress:
    texture = []
    microstructure = []
    otherPoint = -1                                                                                 # ensure to create first microstructure
    matPoints = np.zeros(info['grid'].prod(),dtype='i')                                             # index of microstructure in geom file
    for myPoint in xrange(info['grid'].prod()):
      myTexture = -1
      for otherPoint in xrange(len(microstructure)):
        otherEulers = eulerangles[texture[microstructure[otherPoint][0]]]
        otherPhase  = microstructure[otherPoint][1]
        if all(abs(eulerangles[myPoint]-otherEulers)<1e-6) and phase[myPoint] == otherPhase:        # common microstructure
           matPoints[myPoint] = otherPoint+1                                                        # use other point's microstructure, +1 because starting with 1 (.config) instead of 0 (python)
           otherPoint = -2                                                                          # never create new microstructure
           break
        elif all(eulerangles[myPoint] == otherEulers):                                              # found common texture and store it
           myTexture = microstructure[otherPoint][0]
      if otherPoint == len(microstructure)-1:                                                       # did not found matching microstructure
        if myTexture == -1:                                                                         # did not even found matching texture
          myTexture = len(texture)
          texture.append(myPoint)
        microstructure.append([myTexture,phase[myPoint]])
        matPoints[myPoint] = len(microstructure)                                                    # use the new microstructure
  else:
    texture = [i for i in xrange(info['grid'][0]*info['grid'][1])]
    microstructure = [[i,phase[i]] for i in xrange(info['grid'][0]*info['grid'][1])]
    
  formatOut = 1+int(math.log10(len(texture)))
  textureOut =['\n\n<texture>']
  for i in xrange(len(texture)):
    textureOut +=       ['[Texture%s]\n'%str(texture[i]+1).zfill(formatOut) + \
                          'axes %s %s %s\n'%(options.axes[0],options.axes[1],options.axes[2]) +\
                          '(gauss)\tphi1 %4.2f\tPhi %4.2f\tphi2 %4.2f\tscatter 0.0\tfraction 1.0\n'%tuple(eulerangles[texture[i],...])
                         ]
  formatOut = 1+int(math.log10(len(microstructure)))
  microstructureOut =['<microstructure>']
  for i in xrange(len(microstructure)):
    microstructureOut += ['[Grain%s]\n'%str(i+1).zfill(formatOut) + \
                         'crystallite\t%i\n'%options.crystallite + \
                         '(constituent)\tphase %i\ttexture %i\tfraction 1.0\n'%(microstructure[i][1],microstructure[i][0]+1)
                         ]

  info['microstructures'] = len(microstructure)
  info['size'] = step[0]*info['grid'][0],step[1]*info['grid'][1],min(step)

#--- report ---------------------------------------------------------------------------------------
  file['croak'].write('grid     a b c:  %s\n'%(' x '.join(map(str,info['grid']))) + \
                      'size     x y z:  %s\n'%(' x '.join(map(str,info['size']))) + \
                      'origin   x y z:  %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization:  %i\n'%info['homogenization'] + \
                      'microstructures: %i\n\n'%info['microstructures'])

  if np.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    continue
  if np.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    continue


#--- write data -----------------------------------------------------------------------------------
  if options.config:
    file['output'].write('\n'.join(microstructureOut+ textureOut) + '\n')
  else:
    header = [' '.join([scriptID] + sys.argv[1:]),
              "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
              "size\tx %f\ty %f\tz %f"%(info['size'][0],info['size'][1],info['size'][2],),
              "origin\tx %f\ty %f\tz %f"%(info['origin'][0],info['origin'][1],info['origin'][2],),
              "microstructures\t%i"%info['microstructures'],
              "homogenization\t%i"%info['homogenization'],
              ]
    file['output'].write('\n'.join(['%i\theader'%(len(header))] + header) + '\n')
    if options.compress:
      matPoints = matPoints.reshape((info['grid'][1],info['grid'][0]))
      np.savetxt(file['output'],matPoints,fmt='%0'+str(1+int(math.log10(np.amax(matPoints))))+'d')
    else:
      file['output'].write("1 to %i\n"%(info['microstructures']))
  
#--- output finalization -------------------------------------------------------------------------- 
  if file['name'] != 'STDIN':
    file['output'].close()
    os.rename(file['name']+'_tmp',
              os.path.splitext(file['name'])[0] +'%s'%('_material.config' if options.config else '.geom'))
