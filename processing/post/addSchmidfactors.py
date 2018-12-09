#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

slipSystems = {
'fcc':
    np.array([
    # Slip direction     Plane normal
     [ 0, 1,-1,     1, 1, 1, ],
     [-1, 0, 1,     1, 1, 1, ],
     [ 1,-1, 0,     1, 1, 1, ],
     [ 0,-1,-1,    -1,-1, 1, ],
     [ 1, 0, 1,    -1,-1, 1, ],
     [-1, 1, 0,    -1,-1, 1, ],
     [ 0,-1, 1,     1,-1,-1, ],
     [-1, 0,-1,     1,-1,-1, ],
     [ 1, 1, 0,     1,-1,-1, ],
     [ 0, 1, 1,    -1, 1,-1, ],
     [ 1, 0,-1,    -1, 1,-1, ],
     [-1,-1, 0,    -1, 1,-1, ],
    ],'f'),
'bcc':
    np.array([
    # Slip system <111>{110} 
     [ 1,-1, 1,     0, 1, 1, ],
     [-1,-1, 1,     0, 1, 1, ],
     [ 1, 1, 1,     0,-1, 1, ],
     [-1, 1, 1,     0,-1, 1, ],
     [-1, 1, 1,     1, 0, 1, ],
     [-1,-1, 1,     1, 0, 1, ],
     [ 1, 1, 1,    -1, 0, 1, ],
     [ 1,-1, 1,    -1, 0, 1, ],
     [-1, 1, 1,     1, 1, 0, ],
     [-1, 1,-1,     1, 1, 0, ],
     [ 1, 1, 1,    -1, 1, 0, ],
     [ 1, 1,-1,    -1, 1, 0, ],
    # Slip system <111>{112}
     [-1, 1, 1,     2, 1, 1, ],
     [ 1, 1, 1,    -2, 1, 1, ],
     [ 1, 1,-1,     2,-1, 1, ],
     [ 1,-1, 1,     2, 1,-1, ],
     [ 1,-1, 1,     1, 2, 1, ],
     [ 1, 1,-1,    -1, 2, 1, ],
     [ 1, 1, 1,     1,-2, 1, ],
     [-1, 1, 1,     1, 2,-1, ],
     [ 1, 1,-1,     1, 1, 2, ],
     [ 1,-1, 1,    -1, 1, 2, ],
     [-1, 1, 1,     1,-1, 2, ],
     [ 1, 1, 1,     1, 1,-2, ],
    ],'f'),
'hex':
    np.array([
    # Basal systems <11.0>{00.1} (independent of c/a-ratio, Bravais notation (4 coordinate base))
     [ 2, -1, -1,  0,     0,  0,  0,  1, ],
     [-1,  2, -1,  0,     0,  0,  0,  1, ],
     [-1, -1,  2,  0,     0,  0,  0,  1, ],
    # 1st type prismatic systems <11.0>{10.0}  (independent of c/a-ratio)
     [ 2, -1, -1,  0,     0,  1, -1,  0, ],
     [-1,  2, -1,  0,    -1,  0,  1,  0, ],
     [-1, -1,  2,  0,     1, -1,  0,  0, ],
    # 2nd type prismatic systems <10.0>{11.0} -- a slip; plane normals independent of c/a-ratio
     [ 0,  1,  -1, 0,     2, -1, -1,  0, ],
     [-1,  0,  1,  0,    -1,  2, -1,  0, ],
     [ 1, -1,  0,  0,    -1, -1,  2,  0, ],
    # 1st type 1st order pyramidal systems <11.0>{-11.1} -- plane normals depend on the c/a-ratio
     [ 2, -1, -1,  0,     0,  1, -1,  1, ],
     [-1,  2, -1,  0,    -1,  0,  1,  1, ],
     [-1, -1,  2,  0,     1, -1,  0,  1, ],
     [ 1,  1, -2,  0,    -1,  1,  0,  1, ],
     [-2,  1,  1,  0,     0, -1,  1,  1, ],
     [ 1, -2,  1,  0,     1,  0, -1,  1, ],
    # pyramidal system: c+a slip <11.3>{-10.1} -- plane normals depend on the c/a-ratio
     [ 2, -1, -1,  3,    -1,  1,  0,  1, ],
     [ 1, -2,  1,  3,    -1,  1,  0,  1, ],
     [-1, -1,  2,  3,     1,  0, -1,  1, ],
     [-2,  1,  1,  3,     1,  0, -1,  1, ],
     [-1,  2, -1,  3,     0, -1,  1,  1, ],
     [ 1,  1, -2,  3,     0, -1,  1,  1, ],
     [-2,  1,  1,  3,     1, -1,  0,  1, ],
     [-1,  2, -1,  3,     1, -1,  0,  1, ],
     [ 1,  1, -2,  3,    -1,  0,  1,  1, ],
     [ 2, -1, -1,  3,    -1,  0,  1,  1, ],
     [ 1, -2,  1,  3,     0,  1, -1,  1, ],
     [-1, -1,  2,  3,     0,  1, -1,  1, ],
    # pyramidal system: c+a slip <11.3>{-1-1.2} -- as for hexagonal ice (Castelnau et al. 1996, similar to twin system found below) 
     [ 2, -1, -1,  3,    -2,  1,  1,  2, ], # sorted according to similar twin system
     [-1,  2, -1,  3,     1, -2,  1,  2, ], # <11.3>{-1-1.2} shear = 2((c/a)^2-2)/(3 c/a)
     [-1, -1,  2,  3,     1,  1, -2,  2, ],
     [-2,  1,  1,  3,     2, -1, -1,  2, ],
     [ 1, -2,  1,  3,    -1,  2, -1,  2, ],
     [ 1,  1, -2,  3,    -1, -1,  2,  2, ],
     ],'f'),
}

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add columns listing Schmid factors (and optional trace vector of selected system) for given Euler angles.

""", version = scriptID)

latticeChoices = ('fcc','bcc','hex')
parser.add_option('-l','--lattice',
                  dest = 'lattice', type = 'choice', choices = latticeChoices, metavar='string',
                  help = 'type of lattice structure [%default] {}'.format(latticeChoices))
parser.add_option('--covera',
                  dest = 'CoverA', type = 'float', metavar = 'float',
                  help = 'C over A ratio for hexagonal systems')
parser.add_option('-f', '--force',
                  dest = 'force',
                  type = 'float', nargs = 3, metavar = 'float float float',
                  help = 'force direction in lab frame [%default]')
parser.add_option('-n', '--normal',
                  dest = 'normal',
                  type = 'float', nargs = 3, metavar = 'float float float',
                  help = 'stress plane normal in lab frame [%default]')
parser.add_option('-q', '--quaternion',
                  dest = 'quaternion',
                  metavar = 'string',
                  help = 'quaternion label')

parser.set_defaults(force = (0.0,0.0,1.0),
                    normal = None,
                    lattice = latticeChoices[0],
                    CoverA = math.sqrt(8./3.),
                   )

(options, filenames) = parser.parse_args()

force = np.array(options.force)
force /= np.linalg.norm(force)

if options.normal:
  normal = np.array(options.normal)
  normal /= np.linalg.norm(normal)
  if abs(np.dot(force,normal)) > 1e-3:
    parser.error('stress plane normal not orthogonal to force direction')
else:
  normal = force

slip_direction = np.zeros((len(slipSystems[options.lattice]),3),'f')
slip_normal    = np.zeros_like(slip_direction)


if options.lattice in latticeChoices[:2]:
  slip_direction = slipSystems[options.lattice][:,:3]
  slip_normal    = slipSystems[options.lattice][:,3:]
elif options.lattice == latticeChoices[2]:
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

slip_direction /= np.tile(np.linalg.norm(slip_direction,axis=1),(3,1)).T
slip_normal    /= np.tile(np.linalg.norm(slip_normal   ,axis=1),(3,1)).T

# --- loop over input files ------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------
  if not table.label_dimension(options.quaternion) == 4:
    damask.util.croak('input {} does not have dimension 4.'.format(options.quaternion))
    table.close(dismiss = True)                                                                     # close ASCIItable and remove empty file
    continue

  column = table.label_index(options.quaternion)

# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append(['{id}_'
                       'S[{direction[0]:.1g}_{direction[1]:.1g}_{direction[2]:.1g}]'
                       '({normal[0]:.1g}_{normal[1]:.1g}_{normal[2]:.1g})'\
                       .format(       id = i+1,
                                  normal = theNormal,
                               direction = theDirection,
                              ) for i,(theNormal,theDirection) in enumerate(zip(slip_normal,slip_direction))])
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    o = damask.Orientation(quaternion = np.array(list(map(float,table.data[column:column+4]))))

    table.data_append(  np.abs(  np.sum(slip_direction * (o.quaternion * force) ,axis=1) \
                               * np.sum(slip_normal    * (o.quaternion * normal),axis=1)))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
