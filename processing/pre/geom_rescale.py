#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import numpy as np

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

  grid = geom.get_grid()
  size = geom.get_size()

  new_grid = grid if options.grid is None else \
             np.array([int(o*float(n.lower().replace('x',''))) if n.lower().endswith('x') \
                  else int(n) for o,n in zip(grid,options.grid)],dtype=int)

  new_size = size if options.size is None else \
             np.array([o*float(n.lower().replace('x','')) if n.lower().endswith('x') \
                  else float(n) for o,n in zip(size,options.size)],dtype=float)

  geom.scale(new_grid)
  damask.util.croak(geom.update(microstructure = None,size = new_size))
  geom.add_comments(scriptID + ' ' + ' '.join(sys.argv[1:]))
  geom.to_file(sys.stdout if name is None else name,pack=False)
