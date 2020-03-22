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

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [seedfile(s)]', description = """
Writes vtk file for visualization.

""", version = scriptID)

(options, filenames) = parser.parse_args()
if filenames == []: filenames = [None]

for name in filenames:
    damask.util.report(scriptName,name)

    seeds = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
    v = damask.VTK.from_polyData(seeds.get('pos'))
    for label in seeds.shapes.keys():
        if label == 'pos': pass
        v.add(seeds.get(label),label)

    if name:
        v.write(os.path.splitext(name)[0])
    else:
        sys.stdout.write(v.__repr__())
