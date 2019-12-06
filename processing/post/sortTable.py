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
Sort rows by given (or all) column label(s).

Examples:
With coordinates in columns "x", "y", and "z"; sorting with x slowest and z fastest varying index: --label x,y,z.
""", version = scriptID)


parser.add_option('-l','--label',
                  dest   = 'labels',
                  action = 'extend', metavar = '<string LIST>',
                  help   = 'list of column labels (a,b,c,...)')
parser.add_option('-r','--reverse',
                  dest   = 'reverse',
                  action = 'store_true',
                  help   = 'sort in reverse')

parser.set_defaults(reverse = False,
                   )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

if options.labels is None:
    parser.error('no labels specified.')
for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
    table.sort_by(options.labels,not options.reverse)

    table.to_ASCII(sys.stdout if name is None else name)
