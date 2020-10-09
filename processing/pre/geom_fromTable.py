#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser

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

parser.set_defaults(pos= 'pos')

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

for name in filenames:
    damask.util.report(scriptName,name)

    labels = []
    for l in [options.quaternion,options.phase,options.microstructure]:
        if l is not None: labels.append(l)

    t = damask.Table.load(name)
    geom = damask.Geom.from_table(t,options.pos,labels)
    damask.util.croak(geom)

    geom.save_ASCII(sys.stdout if name is None else os.path.splitext(name)[0]+'.geom',compress=False)
