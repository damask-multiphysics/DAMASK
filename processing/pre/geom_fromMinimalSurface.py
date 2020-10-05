#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


minimal_surfaces = list(damask.Geom._minimal_surface.keys())

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [geomfile]', description = """
Generate a bicontinuous structure of given type.

""", version = scriptID)


parser.add_option('-t','--type',
                  dest = 'type',
                  choices = minimal_surfaces, metavar = 'string',
                  help = 'type of minimal surface [primitive] {%s}' %(','.join(minimal_surfaces)))
parser.add_option('-f','--threshold',
                  dest = 'threshold',
                  type = 'float', metavar = 'float',
                  help = 'threshold value defining minimal surface [%default]')
parser.add_option('-g', '--grid',
                  dest = 'grid',
                  type = 'int', nargs = 3, metavar = 'int int int',
                  help = 'a,b,c grid of hexahedral box [%default]')
parser.add_option('-s', '--size',
                  dest = 'size',
                  type = 'float', nargs = 3, metavar = 'float float float',
                  help = 'x,y,z size of hexahedral box [%default]')
parser.add_option('-p', '--periods',
                  dest = 'periods',
                  type = 'int', metavar = 'int',
                  help = 'number of repetitions of unit cell [%default]')
parser.add_option('--m',
                  dest = 'microstructure',
                  type = 'int', nargs = 2, metavar = 'int int',
                  help = 'two microstructure indices to be used [%default]')

parser.set_defaults(type = minimal_surfaces[0],
                    threshold = 0.0,
                    periods = 1,
                    grid = (16,16,16),
                    size = (1.0,1.0,1.0),
                    microstructure = (1,2),
                   )

(options,filename) = parser.parse_args()


name = None if filename == [] else filename[0]
damask.util.report(scriptName,name)

geom=damask.Geom.from_minimal_surface(options.grid,options.size,options.type,options.threshold,
                                      options.periods,options.microstructure)
damask.util.croak(geom)

geom.save_ASCII(sys.stdout if name is None else name,compress=False)
