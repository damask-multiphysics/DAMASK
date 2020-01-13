#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [geomfile(s)]', description = """
Translate geom description into ASCIItable containing position and microstructure.

""", version = scriptID)

(options, filenames) = parser.parse_args()
if filenames == []: filenames = [None]

for name in filenames:
    damask.util.report(scriptName,name)

    geom = damask.Geom.from_file(StringIO(''.join(sys.stdin.read())) if name is None else name)
    damask.util.croak(geom)

    coord0 = damask.grid_filters.cell_coord0(geom.grid,geom.size,geom.origin).reshape((-1,3),order='F')

    comments = geom.comments \
             + [scriptID + ' ' + ' '.join(sys.argv[1:]),
                "grid\ta {}\tb {}\tc {}".format(*geom.grid),
                "size\tx {}\ty {}\tz {}".format(*geom.size),
                "origin\tx {}\ty {}\tz {}".format(*geom.origin),
                "homogenization\t{}".format(geom.homogenization)]

    table = damask.Table(coord0,{'pos':(3,)},comments)
    table.add('microstructure',geom.microstructure.reshape((-1,1)))

    table.to_ASCII(sys.stdout if name is None else \
                   os.path.splitext(name)[0]+'.txt')
