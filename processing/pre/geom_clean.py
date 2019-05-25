#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os
import sys
import numpy as np
import damask
from io import StringIO
from scipy import ndimage
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

def mostFrequent(arr):
  unique, inverse = np.unique(arr, return_inverse=True)
  return unique[np.argmax(np.bincount(inverse))]


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [geomfile(s)]', description = """
Smooth geometry by selecting most frequent microstructure index within given stencil at each location.

""", version=scriptID)


parser.add_option('-s','--stencil',
                  dest = 'stencil',
                  type = 'int', metavar = 'int',
                  help = 'size of smoothing stencil [%default]')

parser.set_defaults(stencil = 3)

(options, filenames) = parser.parse_args()


# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)
  
  if name is None:
    virt_file = StringIO(''.join(sys.stdin.read()))
    geom = damask.Geom.from_file(virt_file)
  else:
    geom = damask.Geom.from_file(name)
  microstructure = geom.microstructure

  microstructure = ndimage.filters.generic_filter(microstructure,mostFrequent,
                                                  size=(options.stencil,)*3).astype(microstructure.dtype)

  geom.microstructure = microstructure
  geom.add_comment(scriptID + ' ' + ' '.join(sys.argv[1:]))
  
  damask.util.croak(geom)
  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
