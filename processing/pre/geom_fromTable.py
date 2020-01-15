#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser

import numpy as np

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Converts ASCII table. Input can be microstructure or orientation (as quaternion). For the latter,
phase information can be given additionally.

""", version = scriptID)

parser.add_option('--coordinates',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'coordinates label (%default)')
parser.add_option('--phase',
                  dest = 'phase',
                  type = 'string', metavar = 'string',
                  help = 'phase label')
parser.add_option('--microstructure',
                  dest = 'microstructure',
                  type = 'string', metavar = 'string',
                  help = 'microstructure label')
parser.add_option('-q', '--quaternion',
                  dest = 'quaternion',
                  type = 'string', metavar='string',
                  help = 'quaternion label')
parser.add_option('--axes',
                  dest = 'axes',
                  type = 'string', nargs = 3, metavar = ' '.join(['string']*3),
                  help = 'orientation coordinate frame in terms of position coordinate frame [+x +y +z]')
parser.add_option('--homogenization',
                  dest = 'homogenization',
                  type = 'int', metavar = 'int',
                  help = 'homogenization index to be used [%default]')


parser.set_defaults(homogenization = 1,
                    pos            = 'pos',
                   )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

if np.sum([options.quaternion     is not None,
           options.microstructure is not None]) != 1:
  parser.error('need either microstructure or quaternion (and optionally phase) as input.')
if options.microstructure is not None and options.phase is not None:
  parser.error('need either microstructure or phase (and mandatory quaternion) as input.')
if options.axes is not None and not set(options.axes).issubset(set(['x','+x','-x','y','+y','-y','z','+z','-z'])):
  parser.error('invalid axes {} {} {}.'.format(*options.axes))

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
    table.sort_by(['{}_{}'.format(i,options.pos) for i in range(3,0,-1)])                           # x fast, y slow
    grid,size,origin = damask.grid_filters.cell_coord0_gridSizeOrigin(table.get(options.pos))

    config_header  = table.comments

    if   options.microstructure:
        microstructure = table.get(options.microstructure).reshape(grid,order='F')

    elif options.quaternion:
        q     = table.get(options.quaternion)
        phase = table.get(options.phase).astype(int) if options.phase else \
                np.ones((table.data.shape[0],1),dtype=int)

        unique,unique_inverse = np.unique(np.hstack((q,phase)),return_inverse=True,axis=0)
        microstructure = unique_inverse.reshape(grid,order='F') + 1

        config_header = ['<texture>']
        for i,data in enumerate(unique):
            ori = damask.Rotation(data[0:4])
            config_header += ['[Grain{}]'.format(i+1),
                              '(gauss)\tphi1 {:.2f}\tPhi {:.2f}\tphi2 {:.2f}'.format(*ori.asEulers(degrees = True)),
                             ]
            if options.axes is not None: config_header += ['axes\t{} {} {}'.format(*options.axes)]

        config_header += ['<microstructure>']
        for i,data in enumerate(unique):
            config_header += ['[Grain{}]'.format(i+1),
                              '(constituent)\tphase {}\ttexture {}\tfraction 1.0'.format(int(data[4]),i+1),
                             ]

    header = [scriptID + ' ' + ' '.join(sys.argv[1:])]\
           + config_header
    geom = damask.Geom(microstructure,size,origin,
                       homogenization=options.homogenization,comments=header)
    damask.util.croak(geom)

    geom.to_file(sys.stdout if name is None else os.path.splitext(name)[0]+'.geom',pack=False)
