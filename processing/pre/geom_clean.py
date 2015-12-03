#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math
import numpy as np
import damask
from scipy import ndimage
from optparse import OptionParser
from collections import defaultdict

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

def mostFrequent(arr):
  d = defaultdict(int)
  for i in arr: d[i] += 1
  return sorted(d.iteritems(), key=lambda x: x[1], reverse=True)[0][0]                              # return value of most frequent microstructure


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Smooth geometry by selecting most frequent microstructure index within given stencil at each location.

""", version=scriptID)


parser.add_option('-s','--stencil',
                  dest = 'stencil',
                  type = 'int', metavar = 'int',
                  help = 'size of smoothing stencil [%default]')

parser.set_defaults(stencil = 3,
                   )

(options, filenames) = parser.parse_args()


# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False,
                              labeled = False)
  except: continue
  damask.util.report(scriptName,name)

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()
  info,extra_header = table.head_getGeom()

  damask.util.croak(['grid     a b c:  {}'.format(' x '.join(map(str,info['grid']))),
                     'size     x y z:  {}'.format(' x '.join(map(str,info['size']))),
                     'origin   x y z:  {}'.format(' : '.join(map(str,info['origin']))),
                     'homogenization:  {}'.format(info['homogenization']),
                     'microstructures: {}'.format(info['microstructures']),
                    ])

  errors = []
  if np.any(info['grid'] < 1):    errors.append('invalid grid a b c.')
  if np.any(info['size'] <= 0.0): errors.append('invalid size x y z.')
  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --- read data ------------------------------------------------------------------------------------

  microstructure = table.microstructure_read(info['grid']).reshape(info['grid'],order='F')                    # read microstructure

# --- do work ------------------------------------------------------------------------------------

  microstructure = ndimage.filters.generic_filter(microstructure,mostFrequent,size=(options.stencil,options.stencil,options.stencil))

# --- write header ---------------------------------------------------------------------------------

  table.info_clear()
  table.info_append([scriptID + ' ' + ' '.join(sys.argv[1:]),])
  table.head_putGeom(info)
  table.info_append([extra_header])
  table.labels_clear()
  table.head_write()

# --- write microstructure information ------------------------------------------------------------

  formatwidth = int(math.floor(math.log10(microstructure.max())+1))
  table.data = microstructure.reshape((info['grid'][0],np.prod(info['grid'][1:])),order='F').transpose()
  table.data_writeArray('%%%ii'%(formatwidth),delimiter = ' ')

# --- output finalization --------------------------------------------------------------------------

  table.close()                                                                                     # close ASCII table
