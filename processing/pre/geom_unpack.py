#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os
import sys
import math
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Unpack geometry files containing ranges "a to b" and/or "n of x" multiples (exclusively in one line).

""", version = scriptID)

parser.add_option('-1', '--onedimensional',
                  dest   = 'oneD',
                  action = 'store_true',
                  help   = 'output geom file with one-dimensional data arrangement')

parser.set_defaults(oneD = False,
                   )

(options, filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  if name is None:
    virt_file = StringIO(''.join(sys.stdin.read()))
    geom = damask.Geom.from_file(virt_file)
  else:
    geom = damask.Geom.from_file(name)

  damask.util.report(scriptName,name)

  geom.add_comment(scriptID + ' ' + ' '.join(sys.argv[1:]))
  
  damask.util.croak('\n'.join(geom.info()))
  
  if name is None:
    sys.stdout.write(str(geom))
  else:
    geom.to_file(name)
