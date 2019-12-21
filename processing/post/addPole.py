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
Add coordinates of stereographic projection of given direction (pole) in crystal frame.

""", version = scriptID)

parser.add_option('-p',
                  '--pole',
                  dest = 'pole',
                  type = 'float', nargs = 3, metavar = 'float float float',
                  help = 'crystal frame direction for pole figure [%default]')
parser.add_option('--polar',
                  dest = 'polar',
                  action = 'store_true',
                  help = 'output polar coordinates (r,φ) instead of Cartesian coordinates (x,y)')
parser.add_option('-o',
                  '--orientation',
                  dest = 'quaternion',
                  metavar = 'string',
                  help = 'label of crystal orientation given as unit quaternion [%default]')

parser.set_defaults(pole = (1.0,0.0,0.0),
                    quaternion = 'orientation',
                   )

(options, filenames) = parser.parse_args()
if filenames == []: filenames = [None]

pole = np.array(options.pole)
pole /= np.linalg.norm(pole)

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
    orientation = table.get(options.quaternion)
    poles = np.empty((orientation.shape[0],2))
    for i,o in enumerate(orientation):
        rotatedPole = damask.Rotation(o)*pole                                                       # rotate pole according to crystal orientation
        (x,y) = rotatedPole[0:2]/(1.+abs(pole[2]))                                                  # stereographic projection
        poles[i] = [np.sqrt(x*x+y*y),np.arctan2(y,x)] if options.polar else [x,y]                    # cartesian coordinates
  
    table.add('pole_{}{}{}'.format(*options.pole),
              poles,
              scriptID+' '+' '.join(sys.argv[1:]))
    table.to_ASCII(sys.stdout if name is None else name)
