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
Add column(s) containing gradient of requested column(s).
Operates on periodic ordered three-dimensional data sets of scalar and vector fields.
""", version = scriptID)

parser.add_option('-p','--pos','--periodiccellcenter',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'label of coordinates [%default]')
parser.add_option('-l','--label',
                  dest = 'labels',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'label(s) of field values')

parser.set_defaults(pos = 'pos',
                   )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

if options.labels is None: parser.error('no data column specified.')

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
    grid,size,origin = damask.grid_filters.cell_coord0_2_DNA(table.get(options.pos),True)

    for label in options.labels:
        field = table.get(label)
        shape = (1,) if np.prod(field.shape)//np.prod(grid) == 1 else (3,)                          # scalar or vector
        field = field.reshape(np.append(grid[::-1],shape))
        table.add('gradFFT({})'.format(label),
                  damask.grid_filters.gradient(size[::-1],field).reshape((-1,np.prod(shape)*3)),
                  scriptID+' '+' '.join(sys.argv[1:]))
    
    table.to_ASCII(sys.stdout if name is None else name)
