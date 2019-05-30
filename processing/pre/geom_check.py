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
Writes vtk file for visualization.

""", version = scriptID)

(options, filenames) = parser.parse_args()


if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)

  geom = damask.Geom.from_file(StringIO(''.join(sys.stdin.read())) if name is None else name)

  damask.util.croak(geom)
  
  if name is None:
    sys.stdout.write(geom.to_vtk())
  else:
    geom.to_vtk(os.path.splitext(name)[0])
