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
Create regular voxel grid from points in an ASCIItable.

""", version = scriptID)

parser.add_option('-m',
                  '--mode',
                  dest = 'mode',
                  metavar='string',
                  type = 'choice', choices = ['cell','point'],
                  help = 'cell-centered or point-centered coordinates')
parser.add_option('-p',
                  '--pos', '--position',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'label of coordinates [%default]')

parser.set_defaults(mode   = 'cell',
                    pos    = 'pos',
                   )

(options, filenames) = parser.parse_args()
if filenames == []: filenames = [None]

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)

    if   options.mode == 'cell':
        grid, size, origin = damask.grid_filters.cell_coord0_gridSizeOrigin(table.get(options.pos))
    elif options.mode == 'point':
        grid, size, origin = damask.grid_filters.node_coord0_gridSizeOrigin(table.get(options.pos))

    v = damask.VTK.from_rectilinearGrid(grid,size,origin)

    if name:
        v.write('{}_{}({})'.format(os.path.splitext(name)[0],options.pos,options.mode))
    else:
        sys.stdout.write(v.__repr__())
