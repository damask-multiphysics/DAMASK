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
Add column(s) containing eigenvalues and eigenvectors of requested symmetric tensor column(s).

""", version = scriptID)

parser.add_option('-t','--tensor',
                  dest = 'tensor',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'heading of columns containing tensor field values')
parser.add_option('--no-check',
                  dest = 'rh',
                  action = 'store_false',
                  help = 'skip check for right-handed eigenvector basis')

parser.set_defaults(rh = True,
                    )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)

    for tensor in options.tensor:

        t = table.get(tensor).reshape(-1,3,3)
        (u,v) = np.linalg.eigh(damask.mechanics.symmetric(t))
        if options.rh: v[np.linalg.det(v) < 0.0,:,2] *= -1.0

        for i,o in enumerate(['Min','Mid','Max']):
            table.add('eigval{}({})'.format(o,tensor),u[:,i],  scriptID+' '+' '.join(sys.argv[1:]))

        for i,o in enumerate(['Min','Mid','Max']):
            table.add('eigvec{}({})'.format(o,tensor),v[:,:,i],scriptID+' '+' '.join(sys.argv[1:]))
        
    table.to_ASCII(sys.stdout if name is None else name)
