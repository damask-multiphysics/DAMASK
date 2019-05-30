#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser

import numpy as np

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [geomfile(s)]', description = """
Generate homogeneous geometry.

""", version = scriptID)

parser.add_option('-g','--grid',
                  dest = 'grid',
                  type = 'int', nargs = 3, metavar = ' '.join(['int']*3),
                  help = 'a,b,c grid of hexahedral box %default')
parser.add_option('-s', '--size',
                  dest = 'size',
                  type = 'float', nargs = 3, metavar = ' '.join(['float']*3),
                  help = 'x,y,z of geometry size %default')
parser.add_option('-o','--origin',
                  dest = 'origin',
                  type = 'float', nargs = 3, metavar = ' '.join(['float']*3),
                  help = 'x,y,z of geometry origin %default')
parser.add_option('--homogenization',
                  dest = 'homogenization',
                  type = 'int', metavar = 'int',
                  help = 'homogenization index [%default]')
parser.add_option('-f','--fill',
                  dest = 'fill',
                  type = 'float', metavar = 'int',
                  help = 'microstructure index [%default]')

parser.set_defaults(grid = (16,16,16),
                    size = (1.,1.,1.),
                    origin = (0.,0.,0.),
                    homogenization = 1,
                    fill = 1,
                   )

(options, filename) = parser.parse_args()


name = None if filename == [] else filename[0]
damask.util.report(scriptName,name)

dtype = float if np.isnan(options.fill) or int(options.fill) != options.fill else int
geom = damask.Geom(microstructure=np.full(options.grid,options.fill,dtype=dtype),
                   size=options.size,
                   origin=options.origin,
                   homogenization=options.homogenization,
                   comments=scriptID + ' ' + ' '.join(sys.argv[1:]))
damask.util.croak(geom)
  
if name is None:
  sys.stdout.write(str(geom.show()))
else:
  geom.to_file(name)
