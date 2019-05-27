#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import numpy as np

import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [geomfile(s)]', description = """
Increases or decreases the (three-dimensional) canvas.
Grid can be given as absolute or relative values, e.g. 16 16 16 or 2x 0.5x 32.

""", version = scriptID)

parser.add_option('-g','--grid',
                  dest = 'grid',
                  type = 'string', nargs = 3, metavar = ' '.join(['string']*3),
                  help = 'a,b,c grid of hexahedral box')
parser.add_option('-o','--offset',
                  dest = 'offset',
                  type = 'int', nargs = 3, metavar = ' '.join(['int']*3),
                  help = 'a,b,c offset from old to new origin of grid [%default]')
parser.add_option('-f','--fill',
                  dest = 'fill',
                  type = 'float', metavar = 'int',
                  help = 'background microstructure index, defaults to max microstructure index + 1')

parser.set_defaults(offset = (0,0,0))

(options, filenames) = parser.parse_args()


if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)
  
  if name is None:
    virt_file = StringIO(''.join(sys.stdin.read()))
    geom = damask.Geom.from_file(virt_file)
  else:
    geom = damask.Geom.from_file(name)
  microstructure = geom.get_microstructure()
  
  grid = geom.get_grid()
  if options.grid is not None:
    for i,g in enumerate(options.grid):
      grid[i] = int(round(grid[i]*float(g.lower().replace('x','')))) if g.lower().endswith('x') \
           else int(options.grid[i])

  new  = np.full(grid,options.fill if options.fill is not None else np.nanmax(microstructure)+1,
                 microstructure.dtype)
  
  for x in range(microstructure.shape[0]):
    X = x + options.offset[0]
    if not 0 <= X < new.shape[0]: continue 
    for y in range(microstructure.shape[1]):
      Y = y + options.offset[1]
      if not 0 <= Y < new.shape[1]: continue 
      for z in range(microstructure.shape[2]):
        Z = z + options.offset[2]
        if not 0 <= Z < new.shape[2]: continue
        new[X,Y,Z] = microstructure[x,y,z]

  damask.util.croak(geom.update(new,rescale=True))
  geom.add_comment(scriptID + ' ' + ' '.join(sys.argv[1:]))

  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
