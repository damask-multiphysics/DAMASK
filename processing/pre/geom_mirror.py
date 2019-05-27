#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import numpy as np

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [geomfile(s)]', description = """
Mirror along given directions.

""", version=scriptID)

parser.add_option('-d','--direction',
                  dest = 'directions',
                  action = 'extend', metavar = '<string LIST>',
                  help = "directions in which to mirror {'x','y','z'}")
parser.add_option(    '--double',
                  dest = 'double',
                  action = 'store_true',
                  help = 'double the outer layers in mirror direction')             

(options, filenames) = parser.parse_args()

if options.directions is None:
  parser.error('no direction given.')

validDirections = ['x','y','z']  
if not set(options.directions).issubset(validDirections):
  invalidDirections = [str(e) for e in set(options.directions).difference(validDirections)]
  parser.error('invalid directions {}. '.format(*invalidDirections))

limits = [-2,0] if not options.double else [None,None]


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
    microstructure = np.concatenate([microstructure,microstructure[:,:,limits[0]:limits[1]:-1]],2)
  if 'y' in options.directions:
    microstructure = np.concatenate([microstructure,microstructure[:,limits[0]:limits[1]:-1,:]],1)
  if 'x' in options.directions:
    microstructure = np.concatenate([microstructure,microstructure[limits[0]:limits[1]:-1,:,:]],0)
  
  damask.util.croak(geom.update(microstructure,rescale=True))
  geom.add_comment(scriptID + ' ' + ' '.join(sys.argv[1:]))
  
  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
