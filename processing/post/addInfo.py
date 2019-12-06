#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Add info lines to ASCIItable header.

""", version = scriptID)

parser.add_option('-i',
                  '--info',
                  dest = 'info', action = 'extend', metavar = '<string LIST>',
                  help    = 'items to add')

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

if options.info is None:
  parser.error('no info specified.')

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
    table.comments += options.info

    table.to_ASCII(sys.stdout if name is None else name)
