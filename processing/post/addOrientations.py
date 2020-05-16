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
Add quaternion and/or Bunge Euler angle representation of crystal lattice orientation.
Orientation is given by quaternion, Euler angles, rotation matrix, or crystal frame coordinates
(i.e. component vectors of rotation matrix).
Additional (globally fixed) rotations of the lab frame and/or crystal frame can be applied.

""", version = scriptID)

representations = ['quaternion', 'rodrigues', 'eulers', 'matrix', 'axisangle']


parser.add_option('-o',
                  '--output',
                  dest = 'output',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'output orientation formats {{{}}}'.format(', '.join(representations)))
parser.add_option('-d',
                  '--degrees',
                  dest = 'degrees',
                  action = 'store_true',
                  help = 'all angles in degrees')
parser.add_option('-R',
                  '--labrotation',
                  dest='labrotation',
                  type = 'float', nargs = 4, metavar = ' '.join(['float']*4),
                  help = 'axis and angle of additional lab frame rotation [%default]')
parser.add_option('-r',
                  '--crystalrotation',
                  dest='crystalrotation',
                  type = 'float', nargs = 4, metavar = ' '.join(['float']*4),
                  help = 'axis and angle of additional crystal frame rotation [%default]')
parser.add_option('--eulers',
                  dest = 'eulers',
                  metavar = 'string',
                  help = 'Euler angles label')
parser.add_option('--rodrigues',
                  dest = 'rodrigues',
                  metavar = 'string',
                  help = 'Rodrigues vector label')
parser.add_option('--matrix',
                  dest = 'matrix',
                  metavar = 'string',
                  help = 'orientation matrix label')
parser.add_option('--quaternion',
                  dest = 'quaternion',
                  metavar = 'string',
                  help = 'quaternion label')
parser.add_option('-x',
                  dest = 'x',
                  metavar = 'string',
                  help = 'label of lab x vector (expressed in crystal coords)')
parser.add_option('-y',
                  dest = 'y',
                  metavar = 'string',
                  help = 'label of lab y vector (expressed in crystal coords)')
parser.add_option('-z',
                  dest = 'z',
                  metavar = 'string',
                  help = 'label of lab z vector (expressed in crystal coords)')
parser.add_option('--lattice',
                  dest = 'lattice',
                  metavar = 'string',
                  help = 'lattice structure to reduce rotation into fundamental zone')

parser.set_defaults(output = [],
                    labrotation     = (1.,1.,1.,0.),                                                # no rotation about (1,1,1)
                    crystalrotation = (1.,1.,1.,0.),                                                # no rotation about (1,1,1)
                    lattice = None,
                   )

(options, filenames) = parser.parse_args()
if filenames == []: filenames = [None]

if options.output == [] or (not set(options.output).issubset(set(representations))):
    parser.error('output must be chosen from {}.'.format(', '.join(representations)))

input = [options.eulers     is not None,
         options.rodrigues  is not None,
         options.x          is not None and \
         options.y          is not None and \
         options.z          is not None,
         options.matrix     is not None,
         options.quaternion is not None,
        ]

if np.sum(input) != 1: parser.error('needs exactly one input format.')

r = damask.Rotation.from_axis_angle(np.array(options.crystalrotation),options.degrees,normalise=True)
R = damask.Rotation.from_axis_angle(np.array(options.labrotation),options.degrees,normalise=True)

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)

    if   options.eulers is not None:
        label = options.eulers
        print(np.max(table.get(options.eulers),axis=0))
        o = damask.Rotation.from_Eulers(table.get(options.eulers), options.degrees)
    elif options.rodrigues is not None:
        label = options.rodrigues
        o = damask.Rotation.from_Rodrigues(table.get(options.rodrigues))
    elif options.matrix is not None:
        label = options.matrix
        o = damask.Rotation.from_matrix(table.get(options.matrix).reshape(-1,3,3))
    elif options.x is not None:
        label = '<{},{},{}>'.format(options.x,options.y,options.z)
        M = np.block([table.get(options.x),table.get(options.y),table.get(options.z)]).reshape(-1,3,3)
        o = damask.Rotation.from_matrix(M/np.linalg.norm(M,axis=0))
    elif options.quaternion is not None:
        label = options.quaternion
        o = damask.Rotation.from_quaternion(table.get(options.quaternion))

    o = r.broadcast_to(o.shape) @ o @ R.broadcast_to(o.shape)

    #if options.lattice is not None:
    #  o = damask.Orientation(rotation = o,lattice = options.lattice).reduced().rotation


    if 'rodrigues' in options.output:
        table.add('ro({})'.format(label),o.as_rodrigues(),                scriptID+' '+' '.join(sys.argv[1:]))
    if 'eulers' in options.output:
        table.add('eu({})'.format(label),o.as_Eulers(options.degrees),    scriptID+' '+' '.join(sys.argv[1:]))
    if 'quaternion' in options.output:
        table.add('qu({})'.format(label),o.as_quaternion(),               scriptID+' '+' '.join(sys.argv[1:]))
    if 'matrix' in options.output:
        table.add('om({})'.format(label),o.as_matrix(),                   scriptID+' '+' '.join(sys.argv[1:]))
    if 'axisangle' in options.output:
        table.add('om({})'.format(label),o.as_axisangle(options.degrees), scriptID+' '+' '.join(sys.argv[1:]))

    table.to_ASCII(sys.stdout if name is None else name)
