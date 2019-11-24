#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [geomfile(s)]', description = """
Mirror along given directions.

""", version=scriptID)

validDirections = ['x','y','z']

parser.add_option('-d','--direction',
                  dest = 'directions',
                  action = 'extend', metavar = '<string LIST>',
                  help = "directions in which to mirror {{{}}}".format(','.join(validDirections)))
parser.add_option(    '--reflect',
                  dest = 'reflect',
                  action = 'store_true',
                  help = 'reflect (include) outermost layers')

parser.set_defaults(reflect = False)

(options, filenames) = parser.parse_args()

if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)

  geom = damask.Geom.from_file(StringIO(''.join(sys.stdin.read())) if name is None else name)
  damask.util.croak(geom.mirror(options.directions,options.reflect))
  geom.add_comments(scriptID + ' ' + ' '.join(sys.argv[1:]))
  geom.to_file(sys.stdout if name is None else name)
