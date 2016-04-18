#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """

Generate geometry description and material configuration from EBSD data in given square-gridded 'ang' file.
Two phases can be discriminated based on threshold value in a given data column.

""", version = scriptID)

parser.add_option('--column',              dest='column', type='int', metavar = 'int',
                  help='data column to discriminate between both phases [%default]')
parser.add_option('-t','--threshold',      dest='threshold', type='float', metavar = 'float',
                  help='threshold value for phase discrimination [%default]')
parser.add_option('--homogenization',      dest='homogenization', type='int', metavar = 'int',
                  help='homogenization index for <microstructure> configuration [%default]')
parser.add_option('--phase',               dest='phase', type='int', nargs = 2, metavar = 'int int',
                  help='phase indices for <microstructure> configuration %default')
parser.add_option('--crystallite',         dest='crystallite', type='int', metavar = 'int',
                  help='crystallite index for <microstructure> configuration [%default]')
parser.add_option('-c','--compress',            dest='compress', action='store_true',
                  help='lump identical microstructure and texture information [%default]')
parser.add_option('-a', '--axes',         dest='axes', nargs = 3, metavar = 'string string string',
                  help='Euler angle coordinate system for <texture> configuration x,y,z = %default')
parser.add_option('-p', '--precision',    dest='precision', choices=['0','1','2','3'], metavar = 'int',
                  help = 'euler angles decimal places for output format and compressing (0,1,2,3) [2]')
                  
parser.set_defaults(column         = 11,
                    threshold      = 0.5,
                    homogenization = 1,
                    phase          = [1,2],
                    crystallite    = 1,
                    compress       = False,
                    axes           = ['y','x','-z'],
                    precision      = '2',
                 )
(options,filenames) = parser.parse_args()


for i in options.axes:
  if i.lower() not in ['x','+x','-x','y','+y','-y','z','+z','-z']:
    parser.error('invalid axes %s %s %s' %(options.axes[0],options.axes[1],options.axes[2]))

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              outname = os.path.splitext(name)[-2]+'.geom' if name else name,
                              buffered = False, labeled = False)
  except: continue
  damask.util.report(scriptName,name)

  info = {
          'grid':   np.ones(3,'i'),
          'size':   np.zeros(3,'d'),
          'origin': np.zeros(3,'d'),
          'microstructures': 0,
          'homogenization':  options.homogenization
         }
         
  coords = [{},{},{1:True}]
  pos = {'min':[ float("inf"), float("inf")],
         'max':[-float("inf"),-float("inf")]}

  phase =       []
  eulerangles = []

# --------------- read data -----------------------------------------------------------------------
  errors  = []
  while table.data_read():
    words = table.data
    if words[0] == '#':                                                                               # process initial comments/header block
      if len(words) > 2:
        if words[2].lower() == 'hexgrid': 
          errors.append('The file has HexGrid format. Please first convert to SquareGrid...')
          break
    else:
      currPos = words[3:5]
      for i in xrange(2):
        coords[i][currPos[i]] = True
      currPos = map(float,currPos)
      for i in xrange(2):
        pos['min'][i] = min(pos['min'][i],currPos[i])
        pos['max'][i] = max(pos['max'][i],currPos[i])
      eulerangles.append(map(math.degrees,map(float,words[:3])))
      phase.append(options.phase[int(float(words[options.column-1]) > options.threshold)])

  if errors  != []:
    damask.util.croak(errors)
    continue

# --------------- determine size and grid ---------------------------------------------------------
  info['grid'] = np.array(map(len,coords),'i')
  info['size'][0:2] = info['grid'][0:2]/(info['grid'][0:2]-1.0)* \
                                        np.array([pos['max'][0]-pos['min'][0],
                                                  pos['max'][1]-pos['min'][1]],'d')
  info['size'][2]=info['size'][0]/info['grid'][0]
  eulerangles = np.array(eulerangles,dtype='f').reshape(info['grid'].prod(),3)
  phase       = np.array(phase,dtype='i').reshape(info['grid'].prod())

  limits = [360,180,360]
  if any([np.any(eulerangles[:,i]>=limits[i]) for i in [0,1,2]]):
    errors.append('Error: euler angles out of bound. Ang file might contain unidexed poins.')
    for i,angle in enumerate(['phi1','PHI','phi2']):
      for n in np.nditer(np.where(eulerangles[:,i]>=limits[i]),['zerosize_ok']):
        errors.append('%s in line %i (%4.2f %4.2f %4.2f)\n'
                                     %(angle,n,eulerangles[n,0],eulerangles[n,1],eulerangles[n,2]))
  if errors  != []:
    damask.util.croak(errors)
    continue

  eulerangles=np.around(eulerangles,int(options.precision))                                         # round to desired precision
# ensure, that rounded euler angles are not out of bounds (modulo by limits)
  for i,angle in enumerate(['phi1','PHI','phi2']):
    eulerangles[:,i]%=limits[i]

# scale angles by desired precision and convert to int. create unique integer key from three euler angles by
# concatenating the string representation with leading zeros and store as integer and search unique euler angle keys.
# Texture IDs are the indices of the first occurrence, the inverse is used to construct the microstructure
# create a microstructure (texture/phase pair) for each point using unique texture IDs.
# Use longInt (64bit, i8) because the keys might be long
  if options.compress:
    formatString='{0:0>'+str(int(options.precision)+3)+'}'
    euleranglesRadInt = (eulerangles*10**int(options.precision)).astype('int')                      
    eulerKeys = np.array([int(''.join(map(formatString.format,euleranglesRadInt[i,:]))) \
                                                             for i in xrange(info['grid'].prod())]) 
    devNull, texture, eulerKeys_idx = np.unique(eulerKeys, return_index = True, return_inverse=True)
    msFull = np.array([[eulerKeys_idx[i],phase[i]] for i in xrange(info['grid'].prod())],'i8')      
    devNull,msUnique,matPoints = np.unique(msFull.view('c16'),True,True)
    matPoints+=1
    microstructure = np.array([msFull[i] for i in msUnique])                                        # pick only unique microstructures
  else:
    texture = np.arange(info['grid'].prod())
    microstructure = np.hstack( zip(texture,phase) ).reshape(info['grid'].prod(),2)                 # create texture/phase pairs
  formatOut = 1+int(math.log10(len(texture)))
  textureOut =['<texture>']

  eulerFormatOut='%%%i.%if'%(int(options.precision)+4,int(options.precision))
  outStringAngles='(gauss) phi1 '+eulerFormatOut+' Phi '+eulerFormatOut+' phi2 '+eulerFormatOut+' scatter 0.0 fraction 1.0'
  for i in xrange(len(texture)):
    textureOut +=       ['[Texture%s]'%str(i+1).zfill(formatOut),
                          'axes %s %s %s'%(options.axes[0],options.axes[1],options.axes[2]),
                          outStringAngles%tuple(eulerangles[texture[i],...])
                         ]
  formatOut = 1+int(math.log10(len(microstructure)))
  microstructureOut =['<microstructure>']
  for i in xrange(len(microstructure)):
    microstructureOut += ['[Grain%s]'%str(i+1).zfill(formatOut),
                         'crystallite\t%i'%options.crystallite,
                         '(constituent)\tphase %i\ttexture %i\tfraction 1.0'%(microstructure[i,1],microstructure[i,0]+1)
                         ]

  info['microstructures'] = len(microstructure)

#--- report ---------------------------------------------------------------------------------------
  damask.util.croak('grid     a b c:  %s\n'%(' x '.join(map(str,info['grid']))) +
              'size     x y z:  %s\n'%(' x '.join(map(str,info['size']))) +
              'origin   x y z:  %s\n'%(' : '.join(map(str,info['origin']))) +
              'homogenization:  %i\n'%info['homogenization'] +
              'microstructures: %i\n\n'%info['microstructures'])

  if np.any(info['grid'] < 1):
    damask.util.croak('invalid grid a b c.\n')
    continue
  if np.any(info['size'] <= 0.0):
    damask.util.croak('invalid size x y z.\n')
    continue


#--- write data/header --------------------------------------------------------------------------------
  table.info_clear()
  table.info_append([' '.join([scriptID] + sys.argv[1:]),
                     "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
                     "size\tx %f\ty %f\tz %f"%(info['size'][0],info['size'][1],info['size'][2],),
                     "origin\tx %f\ty %f\tz %f"%(info['origin'][0],info['origin'][1],info['origin'][2],),
                     "microstructures\t%i"%info['microstructures'],
                     "homogenization\t%i"%info['homogenization'],
                     ] +
                     [line for line in microstructureOut + textureOut]
                    )
  table.head_write()
  if options.compress:
    matPoints = matPoints.reshape((info['grid'][1],info['grid'][0]))
    table.data = matPoints
    table.data_writeArray('%%%ii'%(1+int(math.log10(np.amax(matPoints)))),delimiter=' ')
  else:
    table.output_write("1 to %i\n"%(info['microstructures']))
  table.output_flush()

# --- output finalization --------------------------------------------------------------------------

  table.close()
