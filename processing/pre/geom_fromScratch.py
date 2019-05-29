#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os
import sys
import numpy as np
import damask

from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog option(s) [geomfile(s)]', description = """
Generate homogeneous geometry.

""", version = scriptID)

parser.add_option('-g',
                  '--grid',
                  dest = 'grid',
                  type = 'string', nargs = 3, metavar = ' '.join(['string']*3),
                  help = 'a,b,c grid of hexahedral box %default')
parser.add_option('-s',
                  '--size',
                  dest = 'size',
                  type = 'float', nargs = 3, metavar = ' '.join(['float']*3),
                  help = 'x,y,z of geometry size %default')
parser.add_option('-o',
                  '--origin',
                  dest = 'origin',
                  type = 'float', nargs = 3, metavar = ' '.join(['float']*3),
                  help = 'x,y,z of geometry origin %default')
parser.add_option('--homogenization',
                  dest = 'homogenization',
                  type = 'int', metavar = 'int',
                  help = 'homogenization index [%default]')
parser.add_option('-f',
                  '--fill',
                  dest = 'fill',
                  type = 'float', metavar = 'float',
                  help = 'microstructure index [%default]')

parser.set_defaults(grid = [16,16,16],
                    size = [1.,1.,1.],
                    origin = [0.,0.,0.],
                    homogenization = 1,
                    fill = 1,
                   )

(options, filenames) = parser.parse_args()

dtype = float if np.isnan(options.fill) or int(options.fill) != options.fill else int
grid   = list(map(  int,options.grid))
size   = list(map(float,options.size))
origin = list(map(float,options.origin))
homogenization = int(options.homogenization)
fill = float(options.fill)

damask.util.report(scriptName,'')

geom = damask.Geom(microstructure=(fill * np.ones(grid)).astype(dtype),
                   size=size,
                   origin=origin,
                   homogenization=homogenization,
                   comments=scriptID + ' ' + ' '.join(sys.argv[1:]))

damask.util.croak(geom)
  
sys.stdout.write(str(geom.show()))
