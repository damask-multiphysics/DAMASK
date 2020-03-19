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

parser = OptionParser(option_class=damask.extendableOption,
                      usage='%prog options [ASCIItable(s)]',
                      description = 'Add scalars/vectors, tensors, and/or a RGB tuples from ASCIItable '
                                  + 'to existing VTK file (.vtr/.vtu/.vtp).',
                      version = scriptID)

parser.add_option(      '--vtk',
                  dest = 'vtk',
                  type = 'string', metavar = 'string',
                  help = 'VTK file name')
parser.add_option('-d', '--data',
                  dest = 'data',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'scalar/vector value(s) label(s)')
parser.add_option('-t', '--tensor',
                  dest = 'tensor',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'tensor (3x3) value label(s)')
parser.add_option('-c', '--color',
                  dest = 'color',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'RGB color tuple label')

parser.set_defaults(data = [],
                    tensor = [],
                    color = [],
)

(options, filenames) = parser.parse_args()
if filenames == []: filenames = [None]

if not options.vtk:
    parser.error('No VTK file specified.')

for name in filenames:
    damask.util.report(scriptName,name)
  
    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
    vtk = damask.VTK.from_file(options.vtk)

    for data in options.data+options.tensor:
        vtk.add(table.get(data),data)
    for color in options.color:
        vtk.add((table.get(color)*255).astype(np.uint8),color)
  
    vtk.write(options.vtk)
