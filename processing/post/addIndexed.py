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

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Add data in column(s) of mapped ASCIItable selected from the row indexed by the value in a mapping column.
Row numbers start at 1.

""", version = scriptID)

parser.add_option('--index',
                  dest = 'index',
                  type = 'string', metavar = 'string',
                  help = 'column label containing row index')
parser.add_option('-o','--offset',
                  dest = 'offset',
                  type = 'int', metavar = 'int',
                  help = 'constant offset for index column value [%default]')
parser.add_option('-l','--label',
                  dest = 'label',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'column label(s) to be appended')
parser.add_option('-a','--asciitable',
                  dest = 'asciitable',
                  type = 'string', metavar = 'string',
                  help = 'indexed ASCIItable')

parser.set_defaults(offset = 0,
                   )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

if options.label is None:
  parser.error('no data columns specified.')
if options.index is None:
  parser.error('no index column given.')

for name in filenames:
    damask.util.report(scriptName,name)

    table        = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
    indexedTable = damask.Table.from_ASCII(options.asciitable)
    idx = np.reshape(table.get(options.index).astype(int) + options.offset,(-1))-1
    
    for data in options.label:
        table.add(data+'_addIndexed',indexedTable.get(data)[idx],scriptID+' '+' '.join(sys.argv[1:]))

    table.to_ASCII(sys.stdout if name is None else name)
