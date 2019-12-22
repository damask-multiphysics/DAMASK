#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser

import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Append data of ASCIItable(s) column-wise.

""", version = scriptID)

parser.add_option('-a', '--add','--table',
                  dest = 'table',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'tables to add')

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

if options.table is None:
      parser.error('no table specified.')

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)

    for addTable in options.table:
        table2 = damask.Table.from_ASCII(addTable)
        table2.data = table2.data[:table.data.shape[0]]
        table.join(table2)

    table.to_ASCII(sys.stdout if name is None else name)
