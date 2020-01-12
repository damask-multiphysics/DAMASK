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
Add RGB color value corresponding to TSL-OIM scheme for inverse pole figures.

""", version = scriptID)

parser.add_option('-p',
                  '--pole',
                  dest = 'pole',
                  type = 'float', nargs = 3, metavar = 'float float float',
                  help = 'lab frame direction for inverse pole figure [%default]')
parser.add_option('-s',
                  '--symmetry',
                  dest = 'symmetry',
                  type = 'choice', choices = damask.Symmetry.lattices[1:], metavar='string',
                  help = 'crystal symmetry [%default] {{{}}} '.format(', '.join(damask.Symmetry.lattices[1:])))
parser.add_option('-o',
                  '--orientation',
                  dest = 'quaternion',
                  metavar = 'string',
                  help = 'label of crystal orientation given as unit quaternion [%default]')

parser.set_defaults(pole = (0.0,0.0,1.0),
                    quaternion = 'orientation',
                    symmetry = damask.Symmetry.lattices[-1],
                   )

(options, filenames) = parser.parse_args()
if filenames == []: filenames = [None]

# damask.Orientation requires Bravais lattice, but we are only interested in symmetry
symmetry2lattice={'cubic':'fcc','hexagonal':'hex','tetragonal':'bct'}
lattice = symmetry2lattice[options.symmetry]

pole = np.array(options.pole)
pole /= np.linalg.norm(pole)

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
    orientation = table.get(options.quaternion)
    color = np.empty((orientation.shape[0],3))
    for i,o in enumerate(orientation):
        color[i] = damask.Orientation(o,lattice = lattice).IPFcolor(pole)
 
    table.add('IPF_{:g}{:g}{:g}_{sym}'.format(*options.pole,sym = options.symmetry.lower()),
              color,
              scriptID+' '+' '.join(sys.argv[1:]))
    table.to_ASCII(sys.stdout if name is None else name)
