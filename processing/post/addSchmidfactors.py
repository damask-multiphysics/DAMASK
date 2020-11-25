#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import numpy as np

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

slipSystems = {
'fcc': damask.lattice.kinematics['cF']['slip'][:12],
'bcc': damask.lattice.kinematics['cI']['slip'],
'hex': damask.lattice.kinematics['hP']['slip'],
}

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(usage='%prog options [ASCIItable(s)]', description = """
Add columns listing Schmid factors (and optional trace vector of selected system) for given Euler angles.

""", version = scriptID)

lattice_choices = list(slipSystems.keys())
parser.add_option('-l',
                  '--lattice',
                  dest = 'lattice', type = 'choice', choices = lattice_choices, metavar='string',
                  help = 'type of lattice structure [%default] {}'.format(lattice_choices))
parser.add_option('--covera',
                  dest = 'CoverA', type = 'float', metavar = 'float',
                  help = 'C over A ratio for hexagonal systems [%default]')
parser.add_option('-f',
                  '--force',
                  dest = 'force',
                  type = 'float', nargs = 3, metavar = 'float float float',
                  help = 'force direction in lab frame [%default]')
parser.add_option('-n',
                  '--normal',
                  dest = 'normal',
                  type = 'float', nargs = 3, metavar = 'float float float',
                  help = 'stress plane normal in lab frame, per default perpendicular to the force')
parser.add_option('-o',
                  '--orientation',
                  dest = 'quaternion',
                  metavar = 'string',
                  help = 'label of crystal orientation given as unit quaternion [%default]')

parser.set_defaults(force = (0.0,0.0,1.0),
                    quaternion='orientation',
                    normal = None,
                    lattice = lattice_choices[0],
                    CoverA = np.sqrt(8./3.),
                   )

(options, filenames) = parser.parse_args()
if filenames == []: filenames = [None]

force = np.array(options.force)/np.linalg.norm(options.force)

if options.normal is not None:
    normal = np.array(options.normal)/np.linalg.norm(options.ormal)
    if abs(np.dot(force,normal)) > 1e-3:
          parser.error('stress plane normal not orthogonal to force direction')
else:
    normal = force


if   options.lattice in ['bcc','fcc']:
    slip_direction = slipSystems[options.lattice][:,:3]
    slip_normal    = slipSystems[options.lattice][:,3:]
elif options.lattice == 'hex':
    slip_direction = np.zeros((len(slipSystems['hex']),3),'d')
    slip_normal    = np.zeros_like(slip_direction)
    # convert 4 Miller index notation of hex to orthogonal 3 Miller index notation
    for i in range(len(slip_direction)):
        slip_direction[i] = np.array([slipSystems['hex'][i,0]*1.5,
                                     (slipSystems['hex'][i,0] + 2.*slipSystems['hex'][i,1])*0.5*np.sqrt(3),
                                      slipSystems['hex'][i,3]*options.CoverA,
                                     ])
        slip_normal[i]    = np.array([slipSystems['hex'][i,4],
                                     (slipSystems['hex'][i,4] + 2.*slipSystems['hex'][i,5])/np.sqrt(3),
                                      slipSystems['hex'][i,7]/options.CoverA,
                                     ])

slip_direction /= np.linalg.norm(slip_direction,axis=1,keepdims=True)
slip_normal    /= np.linalg.norm(slip_normal,   axis=1,keepdims=True)

labels = ['S[{direction[0]:.1g}_{direction[1]:.1g}_{direction[2]:.1g}]'
           '({normal[0]:.1g}_{normal[1]:.1g}_{normal[2]:.1g})'\
           .format(normal = theNormal, direction = theDirection,
           ) for theNormal,theDirection in zip(slip_normal,slip_direction)]

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.load(StringIO(''.join(sys.stdin.read())) if name is None else name)

    o = damask.Rotation.from_quaternion(table.get(options.quaternion))

    force  = np.broadcast_to(force, o.shape+(3,))
    normal = np.broadcast_to(normal,o.shape+(3,))
    slip_direction = np.broadcast_to(slip_direction,o.shape+slip_direction.shape)
    slip_normal    = np.broadcast_to(slip_normal,   o.shape+slip_normal.shape)
    S = np.abs(np.einsum('ijk,ik->ij',slip_direction,(o@force))*
               np.einsum('ijk,ik->ij',slip_normal,   (o@normal)))

    for i,label in enumerate(labels):
        table = table.add(label,S[:,i],scriptID+' '+' '.join(sys.argv[1:]))

    table.save((sys.stdout if name is None else name))
