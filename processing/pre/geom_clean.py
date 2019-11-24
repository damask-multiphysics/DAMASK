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
Smooth microstructure by selecting most frequent index within given stencil at each location.

""", version=scriptID)

parser.add_option('-s','--stencil',
                  dest = 'stencil',
                  type = 'int', metavar = 'int',
                  help = 'size of smoothing stencil [%default]')

parser.set_defaults(stencil = 3)

(options, filenames) = parser.parse_args()


if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)

  geom = damask.Geom.from_file(StringIO(''.join(sys.stdin.read())) if name is None else name)

  damask.util.croak(geom.clean(options.stencil))

  geom.add_comments(scriptID + ' ' + ' '.join(sys.argv[1:]))

  geom.to_file(sys.stdout if name is None else name)
