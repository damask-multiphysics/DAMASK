#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser

import h5py
import numpy as np

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [DREAM.3Dfile(s)]', description = """
Converts DREAM.3D file. Input can be cell data (direct pointwise takeover) or grain data (individual
grains are segmented). Requires orientation data as quaternion.

""", version = scriptID)

parser.add_option('-b','--basegroup',
                  dest = 'basegroup',
                  metavar = 'string',
                  help = 'name of the group in "DataContainers" containing the pointwise (and, if applicable grain average) data')
parser.add_option('-p','--pointwise',
                  dest = 'pointwise',
                  metavar = 'string',
                  help = 'name of the group in "DataContainers/<basegroup>" containing pointwise data [%default]')
parser.add_option('-a','--average',
                  dest = 'average',
                  metavar = 'string',
                  help = 'name of the group in "DataContainers</basegroup>" containing grain average data. '\
                       + 'Leave empty for pointwise data')
parser.add_option('--phase',
                  dest = 'phase',
                  type = 'string',
                  metavar = 'string',
                  help = 'name of the dataset containing pointwise/average phase IDs [%default]')
parser.add_option('--microstructure',
                  dest = 'microstructure',
                  type = 'string',
                  metavar = 'string',
                  help = 'name of the dataset connecting pointwise and average data [%default]')
parser.add_option('-q', '--quaternion',
                  dest = 'quaternion',
                  type = 'string',
                  metavar='string',
                  help = 'name of the dataset containing pointwise/average orientation as quaternion [%default]')
parser.add_option('--homogenization',
                  dest = 'homogenization',
                  type = 'int', metavar = 'int',
                  help = 'homogenization index to be used [%default]')

parser.set_defaults(pointwise      = 'CellData',
                    quaternion     = 'Quats',
                    phase          = 'Phases',
                    microstructure = 'FeatureIds',
                    homogenization = 1,
                   )

(options, filenames) = parser.parse_args()

if options.basegroup is None:
  parser.error('No base group selected')

rootDir ='DataContainers'


if filenames == []: parser.error('no input file specified.')

for name in filenames:
  damask.util.report(scriptName,name)

  errors = []
  
  inFile = h5py.File(name, 'r')
  group_geom = os.path.join(rootDir,options.basegroup,'_SIMPL_GEOMETRY')
  try:
    size   = inFile[os.path.join(group_geom,'DIMENSIONS')][...] \
           * inFile[os.path.join(group_geom,'SPACING')][...]
    grid   = inFile[os.path.join(group_geom,'DIMENSIONS')][...]
    origin = inFile[os.path.join(group_geom,'ORIGIN')][...]
  except:
    errors.append('Geometry data ({}) not found'.format(group_geom))
    
    
  group_pointwise = os.path.join(rootDir,options.basegroup,options.pointwise)
  if options.average is None:
    label = 'Point'

    dataset = os.path.join(group_pointwise,options.quaternion)
    try:
      quats = np.reshape(inFile[dataset][...],(np.product(grid),4))
      rot   = [damask.Rotation.fromQuaternion(q,True,P=+1) for q in quats]
    except:
      errors.append('Pointwise orientation (quaternion) data ({}) not readable'.format(dataset))
      
    dataset = os.path.join(group_pointwise,options.phase)
    try:
      phase = np.reshape(inFile[dataset][...],(np.product(grid)))
    except:
      errors.append('Pointwise phase data ({}) not readable'.format(dataset))
    
    microstructure = np.arange(1,np.product(grid)+1,dtype=int).reshape(grid,order='F')

    
  else:
    label = 'Grain'
    
    dataset = os.path.join(group_pointwise,options.microstructure)
    try:
      microstructure = np.transpose(inFile[dataset][...].reshape(grid[::-1]),(2,1,0))               # convert from C ordering
    except:
      errors.append('Link between pointwise and grain average data ({}) not readable'.format(dataset))
    
    group_average = os.path.join(rootDir,options.basegroup,options.average)
    
    dataset = os.path.join(group_average,options.quaternion)
    try:
      rot = [damask.Rotation.fromQuaternion(q,True,P=+1) for q in inFile[dataset][...][1:]]         # skip first entry (unindexed)
    except:
      errors.append('Average orientation data ({}) not readable'.format(dataset))
      
    dataset = os.path.join(group_average,options.phase)
    try:
      phase = [i[0] for i in inFile[dataset][...]][1:]                                              # skip first entry (unindexed)
    except:
      errors.append('Average phase data ({}) not readable'.format(dataset))

  if errors != []:
    damask.util.croak(errors)
    continue

  config_header = ['<texture>']
  for i in range(np.nanmax(microstructure)):
    config_header += ['[{}{}]'.format(label,i+1),
                      '(gauss)\tphi1 {:.2f}\tPhi {:.2f}\tphi2 {:.2f}'.format(*rot[i].asEulers(degrees = True)),
                     ]
  config_header += ['<microstructure>']
  for i in range(np.nanmax(microstructure)):
    config_header += ['[{}{}]'.format(label,i+1),
                      'crystallite 1',
                      '(constituent)\tphase {}\ttexture {}\tfraction 1.0'.format(phase[i],i+1),
                      ]
                       
  header = [scriptID + ' ' + ' '.join(sys.argv[1:])]\
         + config_header
  geom = damask.Geom(microstructure,size,origin,
                     homogenization=options.homogenization,comments=header)
  damask.util.croak(geom)

  geom.to_file(os.path.splitext(name)[0]+'.geom')
