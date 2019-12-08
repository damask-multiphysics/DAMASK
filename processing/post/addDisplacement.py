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
Add displacments resulting from deformation gradient field.
Operates on periodic three-dimensional x,y,z-ordered data sets.
Outputs at cell centers or cell nodes (into separate file).

""", version = scriptID)

parser.add_option('-f',
                  '--defgrad',
                  dest    = 'f',
                  metavar = 'string',
                  help    = 'label of deformation gradient [%default]')
parser.add_option('-p',
                  '--pos', '--position',
                  dest    = 'pos',
                  metavar = 'string',
                  help    = 'label of coordinates [%default]')
parser.add_option('--nodal',
                  dest    = 'nodal',
                  action  = 'store_true',
                  help    = 'output nodal (instead of cell-centered) displacements')

parser.set_defaults(f   = 'f',
                    pos = 'pos',
                   )

(options,filenames) = parser.parse_args()

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
    grid,size,origin = damask.grid_filters.cell_coord0_2_DNA(table.get(options.pos),True)
    
    F = table.get(options.f).reshape(np.append(grid[::-1],(3,3)))
    if options.nodal:
        table = damask.Table(damask.grid_filters.node_coord0(grid[::-1],size[::-1]).reshape((-1,3)),
                             {'pos':(3,)})
        table.add('avg({}).{}'.format(options.f,options.pos),
                  damask.grid_filters.node_displacement_avg(size[::-1],F).reshape((-1,3)),
                  scriptID+' '+' '.join(sys.argv[1:]))
        table.add('fluct({}).{}'.format(options.f,options.pos),
                  damask.grid_filters.node_displacement_fluct(size[::-1],F).reshape((-1,3)),
                  scriptID+' '+' '.join(sys.argv[1:]))
        table.to_ASCII(sys.stdout if name is None else os.path.splitext(name)[0]+'_nodal.txt')
    else:
        table.add('avg({}).{}'.format(options.f,options.pos),
                  damask.grid_filters.cell_displacement_avg(size[::-1],F).reshape((-1,3)),
                  scriptID+' '+' '.join(sys.argv[1:]))
        table.add('fluct({}).{}'.format(options.f,options.pos),
                  damask.grid_filters.cell_displacement_fluct(size[::-1],F).reshape((-1,3)),
                  scriptID+' '+' '.join(sys.argv[1:]))
        table.to_ASCII(sys.stdout if name is None else name)
