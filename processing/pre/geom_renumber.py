#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os
import sys
import numpy as np
from optparse import OptionParser
from io import StringIO
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [file[s]]', description = """
renumber sorted microstructure indices to 1,...,N.

""", version=scriptID)

(options, filenames) = parser.parse_args()

# --- loop over input files ----------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)
  
  if name is None:
    virt_file = StringIO(''.join(sys.stdin.read()))
    geom = damask.Geom.from_file(virt_file)
  else:
    geom = damask.Geom.from_file(name)
  microstructure = geom.microstructure

  grainIDs = np.unique(microstructure)
  renumbered = np.copy(microstructure)
  
  for i, oldID in enumerate(grainIDs):
    renumbered = np.where(microstructure == oldID, i+1, renumbered)

  geom.microstructure = renumbered
  geom.add_comment(scriptID + ' ' + ' '.join(sys.argv[1:]))
  
  damask.util.croak(geom)
  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
