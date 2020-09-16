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
Add cumulative (sum of first to current row) values for given label(s).
""", version = scriptID)

parser.add_option('-l','--label',
                  dest='labels',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'columns to cumulate')
parser.add_option('-p','--product',
                  dest='product', action = 'store_true',
                  help = 'product of values instead of sum')

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

if options.labels is None:
  parser.error('no data column(s) specified.')

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
    for label in options.labels:
        table.add('cum_{}({})'.format('prod'     if options.product else 'sum',label),
                  np.cumprod(table.get(label),0) if options.product else np.cumsum(table.get(label),0),
                  scriptID+' '+' '.join(sys.argv[1:]))

    table.to_file(sys.stdout if name is None else name)
