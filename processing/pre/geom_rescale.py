#!/usr/bin/env python3

import os
import sys
import numpy as np

from io import StringIO
from optparse import OptionParser
from scipy import ndimage

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [geomfile(s)]', description = """
Scales independently in x, y, and z direction in terms of grid and/or size.
Either absolute values or relative factors (like "0.25x") can be used.

""", version = scriptID)

parser.add_option('-g', '--grid',
                  dest = 'grid',
                  type = 'string', nargs = 3, metavar = 'string string string',
                  help = 'a,b,c grid of hexahedral box')
parser.add_option('-s', '--size',
                  dest = 'size',
                  type = 'string', nargs = 3, metavar = 'string string string',
                  help = 'x,y,z size of hexahedral box')

(options, filenames) = parser.parse_args()


if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)
    
  geom = damask.Geom.from_file(StringIO(''.join(sys.stdin.read())) if name is None else name)
  microstructure = geom.get_microstructure()
  grid = geom.get_grid()
  size = geom.get_size()

  new_grid = grid if options.grid is None else \
             np.array([int(o*float(n.lower().replace('x',''))) if n.lower().endswith('x') \
                  else int(n) for o,n in zip(grid,options.grid)],dtype=int)

  new_size = size if options.size is None else \
             np.array([o*float(n.lower().replace('x','')) if n.lower().endswith('x') \
                  else float(n) for o,n in zip(size,options.size)],dtype=float)

  if np.any(new_grid != grid):
    microstructure = ndimage.interpolation.zoom(microstructure, new_grid/grid,output=microstructure.dtype,
                                                order=0,mode='nearest', prefilter=False)

  damask.util.croak(geom.update(microstructure,new_size))
  geom.add_comments(scriptID + ' ' + ' '.join(sys.argv[1:]))

  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
