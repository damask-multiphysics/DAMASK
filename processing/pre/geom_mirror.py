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

validDirections = ['x','y','z']

parser.add_option('-d','--direction',
                  dest = 'directions',
                  action = 'extend', metavar = '<string LIST>',
                  help = "directions in which to mirror {{{}}}".format(','.join(validDirections)))
parser.add_option(    '--periodic',
                  dest = 'periodic',
                  action = 'store_true',
                  help = 'omit periodic copies of outermost layers in mirror direction')             

parser.set_defaults(periodic = False)

(options, filenames) = parser.parse_args()

if options.directions is None:
  parser.error('no direction given.')

if not set(options.directions).issubset(validDirections):
  invalidDirections = [str(e) for e in set(options.directions).difference(validDirections)]
  parser.error('invalid directions {}. '.format(*invalidDirections))

limits = [-2,0] if options.periodic else [None,None]


if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)
  
  geom = damask.Geom.from_file(StringIO(''.join(sys.stdin.read())) if name is None else name)
  
  microstructure = geom.get_microstructure()
  if 'z' in options.directions:
    microstructure = np.concatenate([microstructure,microstructure[:,:,limits[0]:limits[1]:-1]],2)
  if 'y' in options.directions:
    microstructure = np.concatenate([microstructure,microstructure[:,limits[0]:limits[1]:-1,:]],1)
  if 'x' in options.directions:
    microstructure = np.concatenate([microstructure,microstructure[limits[0]:limits[1]:-1,:,:]],0)
  
  damask.util.croak(geom.update(microstructure,rescale=True))
  geom.add_comments(scriptID + ' ' + ' '.join(sys.argv[1:]))
  
  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
