#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os,sys,h5py
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [dream3dfile[s]]', description = """
Convert DREAM3D file to geometry file. This can be done from cell data (direct pointwise takeover) or
from grain data (individual grains are segmented). Requires orientation data as quaternion.

""", version = scriptID)

parser.add_option('-b','--basegroup',
                  dest = 'basegroup', metavar = 'string',
                  help = 'name of the group in "DataContainers" that contains all the data')
parser.add_option('-p','--pointwise',
                  dest = 'pointwise', metavar = 'string',
                  help = 'name of the group in "DataContainers/<basegroup>" that contains pointwise data [%default]')
parser.add_option('-a','--average',
                  dest = 'average', metavar = 'string',
                  help = 'name of the group in "DataContainers</basegroup>" that contains grain average data. '\
                       + 'Leave empty for pointwise data')
parser.add_option('--phase',
                  dest = 'phase',
                  type = 'string', metavar = 'string',
                  help = 'name of the dataset containing pointwise/average phase IDs [%default]')
parser.add_option('--microstructure',
                  dest = 'microstructure',
                  type = 'string', metavar = 'string',
                  help = 'name of the dataset connecting pointwise and average data [%default]')
parser.add_option('-q', '--quaternion',
                  dest = 'quaternion',
                  type = 'string', metavar='string',
                  help = 'name of the dataset containing pointwise/average orientation as quaternion [%default]')

parser.set_defaults(pointwise      = 'CellData',
                    quaternion     = 'Quats',
                    phase          = 'Phases',
                    microstructure = 'FeatureIds',
                    crystallite    = 1,
                   )

(options, filenames) = parser.parse_args()

if options.basegroup is None:
  parser.error('No base group selected')

rootDir ='DataContainers'

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: parser.error('no input file specified.')

for name in filenames:
  try:
    table = damask.ASCIItable(outname = os.path.splitext(name)[0]+'.geom',
                              buffered = False, labeled=False,
                              )
  except: continue
  damask.util.report(scriptName,name)

  errors = []
  
  info = {}
  ori  = []
  inFile = h5py.File(name, 'r')
  group_geom = os.path.join(rootDir,options.basegroup,'_SIMPL_GEOMETRY')
  try:
    info['size']   = inFile[os.path.join(group_geom,'DIMENSIONS')][...] \
                   * inFile[os.path.join(group_geom,'SPACING')][...]
    info['grid']   = inFile[os.path.join(group_geom,'DIMENSIONS')][...]
    info['origin'] = inFile[os.path.join(group_geom,'ORIGIN')][...]
  except:
    errors.append('Geometry data ({}) not found'.format(group_geom))
    
    
  group_pointwise = os.path.join(rootDir,options.basegroup,options.pointwise)
  if options.average is None:
    label = 'point'
    N_microstructure = np.product(info['grid'])
    
    dataset = os.path.join(group_pointwise,options.quaternion)
    try:
      quats = np.reshape(inFile[dataset][...],(N_microstructure,3))
    except:
      errors.append('Pointwise orientation data ({}) not found'.format(dataset))
      
    texture = [damask.Rotation.fromQuaternion(q,P=+1) for q in quats]

    dataset = os.path.join(group_pointwise,options.phase)
    try:
      phase = np.reshape(inFile[dataset][...],(N_microstructure))
    except:
      errors.append('Pointwise phase data ({}) not found'.format(dataset))

    
  else:
    label = 'grain'
    
    dataset = os.path.join(group_pointwise,options.microstructure)
    try:
      microstructure = np.reshape(inFile[dataset][...],(np.product(info['grid'])))
      N_microstructure = np.max(microstructure)
    except:
      errors.append('Link between pointwise and grain average data ({}) not found'.format(dataset))

    group_average = os.path.join(rootDir,options.basegroup,options.average)
    
    dataset = os.path.join(group_average,options.quaternion)
    try:
      texture = [damask.Rotation.fromQuaternion(q,P=+1) for q in inFile[dataset][...][1:]]          # skip first entry (unindexed)
    except:
      errors.append('Average orientation data ({}) not found'.format(dataset))
      
    dataset = os.path.join(group_average,options.phase)
    try:
      phase = [i[0] for i in inFile[dataset][...]][1:]                                              # skip first entry (unindexed)
    except:
      errors.append('Average phase data ({}) not found'.format(dataset))

  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue
    
    
  mat = damask.Material()
  mat.verbose = False

  # dummy <homogenization>
  h = damask.config.material.Homogenization()
  mat.add_section('Homogenization','none',h)
  info['homogenization'] = 1
  
  # <crystallite> placeholder (same for all microstructures at the moment)
  c = damask.config.material.Crystallite()
  mat.add_section('Crystallite','tbd',c)
  
  # <phase> placeholders
  for i in range(np.max(phase)):
    p = damask.config.material.Phase()
    mat.add_section('phase','phase{}-tbd'.format(i+1),p)

  # <texture>
  for i,o in enumerate(texture):
    t = damask.config.material.Texture()
    t.add_component('gauss',{'eulers':o.asEulers(degrees=True)})
    mat.add_section(part='texture', section='{}{}'.format(label,i+1),initialData=t)
  
  # <microstructure>
  for i in range(N_microstructure):
    m = damask.config.material.Microstructure()
    mat.add_section('microstructure','{}{}'.format(label,i+1),m)
    mat.add_microstructure('{}{}'.format(label,i+1),
                           {'phase':  'phase{}-tbd'.format(phase[i]),
                            'texture':'{}{}'.format(label,i+1),
                            'crystallite':'tbd',
                            'fraction':1
                            })
    
  table.info_append([
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {}\tb {}\tc {}".format(*info['grid']),
    "size\tx {}\ty {}\tz {}".format(*info['size']),
    "origin\tx {}\ty {}\tz {}".format(*info['origin']),
    "homogenization\t{}".format(info['homogenization']),
    str(mat).split('\n')
    ])
  table.head_write()
  
  if options.average is None:
    table.data = [1, 'to', format(N_microstructure)]
    table.data_write()
  else:
    table.data = microstructure.reshape(info['grid'][1]*info['grid'][2],info['grid'][0])
    table.data_writeArray()
    
    
  table.close()  
