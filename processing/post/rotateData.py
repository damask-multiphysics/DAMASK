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
Rotate vector and/or tensor column data by given angle around given axis.

""", version = scriptID)

parser.add_option('-d', '--data',
                  dest = 'data',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'vector/tensor value(s) label(s)')
parser.add_option('-r', '--rotation',
                  dest = 'rotation',
                  type = 'float', nargs = 4, metavar = ' '.join(['float']*4),
                  help = 'axis and angle to rotate data [%default]')
parser.add_option('--degrees',
                  dest = 'degrees',
                  action = 'store_true',
                  help = 'angles are given in degrees')

parser.set_defaults(rotation = (1.,1.,1.,0),                                                        # no rotation about (1,1,1)
                    degrees = False,
                   )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

if options.data is None:
    parser.error('no data column specified.')

r = damask.Rotation.from_axis_angle(options.rotation,options.degrees,normalise=True)

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)

    for data in options.data:
        d = table.get(data)
        if table.shapes[data] == (9,): d=d.reshape(-1,3,3)
        d = r.broadcast_to(d.shape[0:1]) @ d
        if table.shapes[data] == (9,): d=d.reshape(-1,9)

        table.set(data,d,scriptID+' '+' '.join(sys.argv[1:]))

    table.to_ASCII(sys.stdout if name is None else name)
