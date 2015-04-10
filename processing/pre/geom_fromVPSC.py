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
Generate geometry description and material configuration from input files used by R.A. Lebensohn.

""", version = scriptID)

parser.add_option('--column',          dest='column', type='int', metavar = 'int',
                  help='data column to discriminate between both phases [%default]')
parser.add_option('-t','--threshold',      dest='threshold', type='float', metavar = 'float',
                  help='threshold value for phase discrimination [%default]')
parser.add_option('--homogenization', dest='homogenization', type='int', metavar = 'int',
                  help='homogenization index for <microstructure> configuration [%default]')
parser.add_option('--phase', dest='phase', type='int', nargs = 2, metavar = 'int int',
                  help='phase indices for <microstructure> configuration %default')
parser.add_option('--crystallite', dest='crystallite', type='int', metavar = 'int',
                  help='crystallite index for <microstructure> configuration [%default]')
parser.add_option('-c', '--configuration', dest='config', action='store_true',
                  help='output material configuration [%default]')
parser.add_option('--compress',            dest='compress', action='store_true',
                  help='lump identical microstructure and texture information [%default]')
parser.add_option('-p', '--precision',    dest='precision', choices=['0','1','2','3'], metavar = 'int',
                  help = 'euler angles decimal places for output format and compressing (0,1,2,3) [2]')

parser.set_defaults(column         = 7)
parser.set_defaults(threshold      = 1.0)
parser.set_defaults(homogenization = 1)
parser.set_defaults(phase          = [1,2])
parser.set_defaults(crystallite    = 1)
parser.set_defaults(config         = False)
parser.set_defaults(compress       = False)
parser.set_defaults(precision      = '2')

(options,filenames) = parser.parse_args()

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
          'grid':   np.zeros(3,'i'),
          'size':   np.zeros(3,'d'),
          'origin': np.zeros(3,'d'),
          'microstructures': 0,
          'homogenization':  options.homogenization
         }

  phase =       []
  eulerangles = []
  point          = 0
  for line in file['input']:
    if line.strip():
      words = line.split()
      currPos = map(float,words[3:6])
      for i in xrange(3):
        if currPos[i] > info['grid'][i]: 
          info['size'][i] = currPos[i]
          info['grid'][i]+=1
      eulerangles.append(map(float,words[:3]))
      phase.append(options.phase[int(float(words[options.column-1]) > options.threshold)])
  
  eulerangles = np.array(eulerangles,dtype='f').reshape(info['grid'].prod(),3)
  phase       = np.array(phase,dtype='i').reshape(info['grid'].prod())
  
  limits = [360,180,360]
  if any([np.any(eulerangles[:,i]>=limits[i]) for i in [0,1,2]]):
    file['croak'].write('Error: euler angles out of bound. Ang file might contain unidexed poins.\n')
    for i,angle in enumerate(['phi1','PHI','phi2']):
      for n in np.nditer(np.where(eulerangles[:,i]>=limits[i]),['zerosize_ok']):
        file['croak'].write('%s in line %i (%4.2f %4.2f %4.2f)\n'
                                     %(angle,n,eulerangles[n,0],eulerangles[n,1],eulerangles[n,2]))
    continue
  eulerangles=np.around(eulerangles,int(options.precision))                                         # round to desired precision
  for i,angle in enumerate(['phi1','PHI','phi2']):
    eulerangles[:,i]%=limits[i]                                                                     # ensure, that rounded euler angles are not out of bounds (modulo by limits)

  if options.compress:
    formatString='{0:0>'+str(int(options.precision)+3)+'}'
    euleranglesRadInt = (eulerangles*10**int(options.precision)).astype('int')                      # scale by desired precision and convert to int
    eulerKeys = np.array([int(''.join(map(formatString.format,euleranglesRadInt[i,:]))) \
                                                             for i in xrange(info['grid'].prod())]) # create unique integer key from three euler angles by concatenating the string representation with leading zeros and store as integer
    devNull, texture, eulerKeys_idx = np.unique(eulerKeys, return_index = True, return_inverse=True)# search unique euler angle keys. Texture IDs are the indices of the first occurence, the inverse is used to construct the microstructure
    msFull = np.array([[eulerKeys_idx[i],phase[i]] for i in xrange(info['grid'].prod())],'i8')      # create a microstructure (texture/phase pair) for each point using unique texture IDs. Use longInt (64bit, i8) because the keys might be long
    devNull,msUnique,matPoints = np.unique(msFull.view('c16'),True,True)
    matPoints+=1
    microstructure = np.array([msFull[i] for i in msUnique])                                        # pick only unique microstructures
  else:
    texture = np.arange(info['grid'].prod())
    microstructure = np.hstack( zip(texture,phase) ).reshape(info['grid'].prod(),2)                 # create texture/phase pairs
  formatOut = 1+int(math.log10(len(texture)))
  textureOut =['\n\n<texture>']

  eulerFormatOut='%%%i.%if'%(int(options.precision)+4,int(options.precision))
  outStringAngles='(gauss) phi1 '+eulerFormatOut+' Phi '+eulerFormatOut+' phi2 '+eulerFormatOut+' scatter 0.0 fraction 1.0\n'
  for i in xrange(len(texture)):
    textureOut +=       ['[Texture%s]\n'%str(i+1).zfill(formatOut) +
                          outStringAngles%tuple(eulerangles[texture[i],...])
                         ]
  formatOut = 1+int(math.log10(len(microstructure)))
  microstructureOut =['<microstructure>']
  for i in xrange(len(microstructure)):
    microstructureOut += ['[Grain%s]\n'%str(i+1).zfill(formatOut) +
                         'crystallite\t%i\n'%options.crystallite +
                         '(constituent)\tphase %i\ttexture %i\tfraction 1.0\n'%(microstructure[i,1],microstructure[i,0]+1)
                         ]

  info['microstructures'] = len(microstructure)

#--- report ---------------------------------------------------------------------------------------
  file['croak'].write('grid     a b c:  %s\n'%(' x '.join(map(str,info['grid']))) +
                      'size     x y z:  %s\n'%(' x '.join(map(str,info['size']))) +
                      'origin   x y z:  %s\n'%(' : '.join(map(str,info['origin']))) +
                      'homogenization:  %i\n'%info['homogenization'] +
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
