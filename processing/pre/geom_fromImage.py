#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
from optparse import OptionParser
from PIL import Image
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """

Generate geometry description from (multilayer) images.
Microstructure index is based on gray scale value (1..256).

""", version = scriptID)

parser.add_option('--homogenization',
                  dest = 'homogenization',
                  type = 'int', metavar = 'int',
                  help = 'homogenization index [%default]')
              
parser.set_defaults(homogenization = 1,
                   )

(options,filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              outname = os.path.splitext(name)[0]+'.geom' if name else name,
                              buffered = False, labeled = False)
  except: continue
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name else ''))

# --- read image ------------------------------------------------------------------------------------

  img = Image.open(name).convert(mode = 'L')                                                        # open and convert to grayscale 8bit

  slice = 0
  while True:
    try:
      img.seek(slice)                                                                               # advance to slice
      layer = np.expand_dims(1+np.array(img,dtype = 'uint16'),axis = 0)                             # read image layer
      microstructure = layer if slice == 0 else np.vstack((microstructure,layer))                   # noqa
      slice += 1                                                                                    # advance to next slice
    except EOFError:
      break

# http://docs.scipy.org/doc/scipy/reference/ndimage.html
# http://scipy-lectures.github.io/advanced/image_processing/

  info = {
          'grid':            np.array(microstructure.shape,'i')[::-1],
          'size':            np.array(microstructure.shape,'d')[::-1],
          'origin':          np.zeros(3,'d'),
          'microstructures': len(np.unique(microstructure)),
          'homogenization':  options.homogenization,
         }

# --- report ---------------------------------------------------------------------------------------

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

# --- write header ---------------------------------------------------------------------------------

  table.info_clear()
  table.info_append([
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=info['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=info['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=info['origin']),
    "homogenization\t{homog}".format(homog=info['homogenization']),
    "microstructures\t{microstructures}".format(microstructures=info['microstructures']),
    ])
  table.labels_clear()
  table.head_write()
  table.output_flush()
  
# --- write microstructure information ------------------------------------------------------------

  formatwidth = int(math.floor(math.log10(microstructure.max())+1))
  table.data = microstructure.reshape((info['grid'][1]*info['grid'][2],info['grid'][0]),order='C')
  table.data_writeArray('%%%ii'%(formatwidth),delimiter = ' ')
    
# --- output finalization --------------------------------------------------------------------------

  table.close()                                                                                     # close ASCII table
