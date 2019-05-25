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

validDirections = ['x','y','z']
parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [geomfile(s)]', description = """
Mirrors spectral geometry description along given directions.

""", version=scriptID)

parser.add_option('-d','--direction',
                  dest = 'directions',
                  action = 'extend', metavar = '<string LIST>',
                  help = "directions in which to mirror {'x','y','z'}")

(options, filenames) = parser.parse_args()

if options.directions is None:
  parser.error('no direction given.')
if not set(options.directions).issubset(validDirections):
  invalidDirections = [str(e) for e in set(options.directions).difference(validDirections)]
  parser.error('invalid directions {}. '.format(*invalidDirections))
  

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)
  
  if name is None:
    virt_file = StringIO(''.join(sys.stdin.read()))
    geom = damask.Geom.from_file(virt_file)
  else:
    geom = damask.Geom.from_file(name)
  microstructure = geom.microstructure

  if 'z' in options.directions:
    microstructure = np.concatenate([microstructure,microstructure[:,:,::-1]],2) # better not double edges
    geom.set_size(geom.get_size()*np.array([1,1,2]))
  if 'y' in options.directions:
    microstructure = np.concatenate([microstructure,microstructure[:,::-1,:]],1) # better not double edges
    geom.set_size(geom.get_size()*np.array([1,2,1]))
  if 'x' in options.directions:
    microstructure = np.concatenate([microstructure,microstructure[::-1,:,:]],0) # better not double edges
    geom.set_size(geom.get_size()*np.array([2,1,1]))
  
  geom.microstructure = microstructure
  geom.add_comment(scriptID + ' ' + ' '.join(sys.argv[1:]))
  
  damask.util.croak(geom)
  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
