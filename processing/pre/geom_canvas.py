#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os
import sys
import numpy as np
import damask

from io import StringIO
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog option(s) [geomfile(s)]', description = """
Changes the (three-dimensional) canvas of a spectral geometry description.
Grid can be given as absolute or relative values, e.g. 16 16 16 or 2x 0.5x 32.

""", version = scriptID)

parser.add_option('-g',
                  '--grid',
                  dest = 'grid',
                  type = 'string', nargs = 3, metavar = ' '.join(['string']*3),
                  help = 'a,b,c grid of hexahedral box. [auto]')
parser.add_option('-o',
                  '--offset',
                  dest = 'offset',
                  type = 'int', nargs = 3, metavar = ' '.join(['int']*3),
                  help = 'a,b,c offset from old to new origin of grid [%default]')
parser.add_option('-f',
                  '--fill',
                  dest = 'fill',
                  type = 'float', metavar = 'float',
                  help = '(background) canvas grain index. "0" selects maximum microstructure index + 1 [%default]')
parser.add_option('--blank',
                  dest = 'blank',
                  action = 'store_true',
                  help = 'blank out (optional) input canvas content')

parser.set_defaults(grid = ['0','0','0'],
                    offset = (0,0,0),
                   )

(options, filenames) = parser.parse_args()

options.grid = ['1','1','1'] if options.blank and options.grid == ['0','0','0'] else options.grid
options.fill = 1 if options.blank and options.fill is None else options.fill

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)
  
  if name is None and options.blank:
    grid = np.array(list(map(int,options.grid)))
    geom = damask.Geom(size=grid,microstructure=options.fill*np.ones(grid))
  else:
    geom = damask.Geom.from_file(StringIO(''.join(sys.stdin.read())) if name is None else name)

  grid = geom.get_grid()
  size = geom.get_size()
  origin = geom.get_origin()
  microstructure = geom.get_microstructure()
  fill = np.nanmax(microstructure)+1 if options.fill is None else options.fill
  dtype = float if np.isnan(fill) or int(fill) != fill or microstructure.dtype==np.float else int

  damask.util.croak(geom)
  
# --- do work ------------------------------------------------------------------------------------
 
  new_grid = np.array([int(o*float(n.lower().replace('x',''))) if n.lower().endswith('x') \
                  else int(n) for o,n in zip(grid,options.grid)],dtype=int)
  new_grid = np.where(new_grid > 0, new_grid,grid)

  microstructure_cropped = np.zeros(new_grid,dtype=dtype)
  microstructure_cropped.fill(fill)
  
  if not options.blank:
    xindex = np.arange(max(options.offset[0],0),min(options.offset[0]+new_grid[0],grid[0]))
    yindex = np.arange(max(options.offset[1],0),min(options.offset[1]+new_grid[1],grid[1]))
    zindex = np.arange(max(options.offset[2],0),min(options.offset[2]+new_grid[2],grid[2]))
    translate_x = [i - options.offset[0] for i in xindex]
    translate_y = [i - options.offset[1] for i in yindex]
    translate_z = [i - options.offset[2] for i in zindex]
    if 0 in map(len,[xindex,yindex,zindex,translate_x,translate_y,translate_z]):
      damask.util.croak('invaldid combination of grid and offset.')
      continue
    microstructure_cropped[min(translate_x):max(translate_x)+1,
                           min(translate_y):max(translate_y)+1,
                           min(translate_z):max(translate_z)+1] \
          = microstructure[min(xindex):max(xindex)+1,
                           min(yindex):max(yindex)+1,
                           min(zindex):max(zindex)+1]

  new_size   = size/grid*new_grid  if np.all(grid > 0) else new_grid
  new_origin = origin + (size/grid if np.all(grid > 0) else new_size/new_grid) * options.offset

  geom.set_microstructure(microstructure_cropped)
  geom.set_size(new_size)
  geom.set_origin(new_origin)
  geom.add_comments(scriptID + ' ' + ' '.join(sys.argv[1:]))

  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
