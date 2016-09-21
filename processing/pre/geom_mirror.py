#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
import damask
from scipy import ndimage
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

validDirections = ['x','y','z']
parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [geomfile(s)]', description = """
Mirrors spectral geometry description along given directions.

""", version=scriptID)

parser.add_option('-d','--direction',
                  dest = 'directions',
                  action = 'extend', metavar = '<string LIST>',
                  help = "directions in which to mirror {'x','y','z'}")

(options, filenames) = parser.parse_args()

if options.directions is None:
  parser.error('no direction given.')
if not set(options.directions).issubset(validDirections):
  invalidDirections = [str(e) for e in set(options.directions).difference(validDirections)]
  parser.error('invalid directions {}. '.format(*invalidDirections))

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False, labeled = False)
  except: continue
  damask.util.report(scriptName,name)

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()
  info,extra_header = table.head_getGeom()

  damask.util.croak(['grid     a b c:  %s'%(' x '.join(map(str,info['grid']))),
                     'size     x y z:  %s'%(' x '.join(map(str,info['size']))),
                     'origin   x y z:  %s'%(' : '.join(map(str,info['origin']))),
                     'homogenization:  %i'%info['homogenization'],
                     'microstructures: %i'%info['microstructures'],
                    ])

  errors = []
  if np.any(info['grid'] < 1):    errors.append('invalid grid a b c.')
  if np.any(info['size'] <= 0.0): errors.append('invalid size x y z.')
  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --- read data ------------------------------------------------------------------------------------

  microstructure = table.microstructure_read(info['grid']).reshape(info['grid'],order='F')          # read microstructure

  if 'z' in options.directions:
    microstructure = np.concatenate([microstructure,microstructure[:,:,::-1]],2)
  if 'y' in options.directions:
    microstructure = np.concatenate([microstructure,microstructure[:,::-1,:]],1)
  if 'x' in options.directions:
    microstructure = np.concatenate([microstructure,microstructure[::-1,:,:]],0)

# --- do work ------------------------------------------------------------------------------------

  newInfo = {
             'size':   microstructure.shape*info['size']/info['grid'],
             'grid':   microstructure.shape,
            }


# --- report ---------------------------------------------------------------------------------------

  remarks = []
  if (any(newInfo['grid']            != info['grid'])):
    remarks.append('--> grid    a b c:  %s'%(' x '.join(map(str,newInfo['grid']))))
  if (any(newInfo['size']            != info['size'])):
    remarks.append('--> size    x y z:  %s'%(' x '.join(map(str,newInfo['size']))))
  if remarks != []: damask.util.croak(remarks)

# --- write header ---------------------------------------------------------------------------------

  table.labels_clear()
  table.info_clear()
  table.info_append([
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=newInfo['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=newInfo['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=info['origin']),
    "homogenization\t{homog}".format(homog=info['homogenization']),
    "microstructures\t{microstructures}".format(microstructures=info['microstructures']),
    extra_header
    ])
  table.head_write()

# --- write microstructure information ------------------------------------------------------------

  formatwidth = int(math.floor(math.log10(microstructure.max())+1))
  table.data = microstructure.reshape((newInfo['grid'][0],np.prod(newInfo['grid'][1:])),order='F').transpose()
  table.data_writeArray('%%%ii'%(formatwidth),delimiter = ' ')

# --- output finalization --------------------------------------------------------------------------

  table.close()                                                                                     # close ASCII table
