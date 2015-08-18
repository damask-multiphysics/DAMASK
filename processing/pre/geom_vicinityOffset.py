#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math
import numpy as np
from scipy import ndimage
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Offset microstructure index for points which see a microstructure different from themselves within a given (cubic) vicinity,
i.e. within the region close to a grain/phase boundary.

""", version = scriptID)

parser.add_option('-v', '--vicinity',
                  dest = 'vicinity',
                  type = 'int', metavar = 'int',
                  help = 'voxel distance checked for presence of other microstructure [%default]')
parser.add_option('-m', '--microstructureoffset',
                  dest='offset',
                  type = 'int', metavar = 'int',
                  help = 'offset (positive or negative) for tagged microstructure indices. '+
                         '"0" selects maximum microstructure index [%default]')

parser.set_defaults(vicinity = 1,
                    offset   = 0,
                   )

(options, filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = ['STDIN']

for name in filenames:
  if not (name == 'STDIN' or os.path.exists(name)): continue
  table = damask.ASCIItable(name = name, outname = name+'_tmp',
                            buffered = False, labeled = False)
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name != 'STDIN' else ''))

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()
  info,extra_header = table.head_getGeom()
  
  table.croak(['grid     a b c:  %s'%(' x '.join(map(str,info['grid']))),
               'size     x y z:  %s'%(' x '.join(map(str,info['size']))),
               'origin   x y z:  %s'%(' : '.join(map(str,info['origin']))),
               'homogenization:  %i'%info['homogenization'],
               'microstructures: %i'%info['microstructures'],
              ])

  errors = []
  if np.any(info['grid'] < 1):    errors.append('invalid grid a b c.')
  if np.any(info['size'] <= 0.0): errors.append('invalid size x y z.')
  if errors != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue

# --- read data ------------------------------------------------------------------------------------

  microstructure = table.microstructure_read(info['grid']).reshape(info['grid'],order='F')          # read microstructure

# --- do work ------------------------------------------------------------------------------------

  newInfo = {
             'microstructures': 0,
            }

  if options.offset == 0: options.offset = microstructure.max()

  microstructure = np.where(ndimage.filters.maximum_filter(microstructure,size=1+2*options.vicinity,mode='wrap') ==
                            ndimage.filters.minimum_filter(microstructure,size=1+2*options.vicinity,mode='wrap'),
                            microstructure, microstructure + options.offset)

  newInfo['microstructures'] = microstructure.max()

# --- report ---------------------------------------------------------------------------------------

  remarks = []
  if (    newInfo['microstructures'] != info['microstructures']): remarks.append('--> microstructures: %i'%newInfo['microstructures'])
  if remarks != []: table.croak(remarks)

# --- write header ---------------------------------------------------------------------------------

  table.labels_clear()
  table.info_clear()
  table.info_append(extra_header+[
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=info['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=info['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=info['origin']),
    "homogenization\t{homog}".format(homog=info['homogenization']),
    "microstructures\t{microstructures}".format(microstructures=newInfo['microstructures']),
    ])
  table.head_write()
  table.output_flush()
  
# --- write microstructure information ------------------------------------------------------------

  formatwidth = int(math.floor(math.log10(microstructure.max())+1))
  table.data = microstructure.reshape((info['grid'][0],info['grid'][1]*info['grid'][2]),order='F').transpose()
  table.data_writeArray('%%%ii'%(formatwidth),delimiter = ' ')
    
# --- output finalization --------------------------------------------------------------------------

  table.close()                                                                                     # close ASCII table
  if name != 'STDIN': os.rename(name+'_tmp',name)                                                   # overwrite old one with tmp new
