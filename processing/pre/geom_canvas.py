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

  geom = damask.Geom.from_file(StringIO(''.join(sys.stdin.read())) if name is None else name)
  origin = geom.get_origin()
  size   = geom.get_size()
  old = new = geom.get_grid()
  offset = np.asarray(options.offset)

  if options.grid is not None:
    new = np.maximum(1,
                     np.array([int(o*float(n.lower().replace('x',''))) if n.lower().endswith('x') \
                          else int(n) for o,n in zip(old,options.grid)],dtype=int))

  canvas = np.full(new,options.fill if options.fill is not None
                  else np.nanmax(geom.microstructure)+1,geom.microstructure.dtype)

  l = np.clip( offset,    0,np.minimum(old  +offset,new))
  r = np.clip( offset+old,0,np.minimum(old*2+offset,new))
  L = np.clip(-offset,    0,np.minimum(new  -offset,old))
  R = np.clip(-offset+new,0,np.minimum(new*2-offset,old))
  canvas[l[0]:r[0],l[1]:r[1],l[2]:r[2]] = geom.microstructure[L[0]:R[0],L[1]:R[1],L[2]:R[2]]


  damask.util.croak(geom.update(canvas,origin=origin+offset*size/old,rescale=True))
  geom.add_comments(scriptID + ' ' + ' '.join(sys.argv[1:]))

  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
