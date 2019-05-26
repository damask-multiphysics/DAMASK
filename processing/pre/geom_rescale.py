#!/usr/bin/env python3

import os
import sys
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
    
  if name is None:
    virt_file = StringIO(''.join(sys.stdin.read()))
    geom = damask.Geom.from_file(virt_file)
  else:
    geom = damask.Geom.from_file(name)
  damask.util.croak(geom)
  microstructure = geom.get_microstructure()

  scale = geom.get_grid().astype('float')
  if options.grid is not None:
    for i,g in enumerate(options.grid):
      scale[i] = scale[i]*float(g.lower().replace('x','')) if g.lower().endswith('x') \
            else float(options.grid[i])/scale[i]
  
  size = geom.get_size()
  if options.size is not None:
    for i,s in enumerate(options.size):
      size[i] = size[i]*float(s.lower().replace('x','')) if s.lower().endswith('x') \
           else options.size[i]

  microstructure = ndimage.interpolation.zoom(microstructure, scale, output=microstructure.dtype,
                                              order=0, mode='nearest', prefilter=False)

  damask.util.croak(geom.update(microstructure,size))
  geom.add_comment(scriptID + ' ' + ' '.join(sys.argv[1:]))

  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
