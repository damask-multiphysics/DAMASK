#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os
import sys
import numpy as np
from io import StringIO
from scipy import ndimage
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [geomfile(s)]', description = """
Rotates spectral geometry description.

""", version=scriptID)

parser.add_option('-r', '--rotation',
                  dest='rotation',
                  type = 'float', nargs = 4, metavar = ' '.join(['float']*4),
                  help = 'rotation given as angle and axis')
parser.add_option('-e', '--eulers',
                  dest = 'eulers',
                  type = 'float', nargs = 3, metavar = ' '.join(['float']*3),
                  help = 'rotation given as Euler angles')
parser.add_option('-d', '--degrees',
                  dest = 'degrees',
                  action = 'store_true',
                  help = 'Euler angles are given in degrees [%default]')
parser.add_option('-m', '--matrix',
                  dest = 'matrix',
                  type = 'float', nargs = 9, metavar = ' '.join(['float']*9),
                  help = 'rotation given as matrix')
parser.add_option('-q', '--quaternion',
                  dest = 'quaternion',
                  type = 'float', nargs = 4, metavar = ' '.join(['float']*4),
                  help = 'rotation given as quaternion')
parser.add_option('-f', '--fill',
                  dest = 'fill',
                  type = 'int', metavar = 'int',
                  help = 'background grain index, defaults to max + 1')

parser.set_defaults(degrees = False,
                   )

(options, filenames) = parser.parse_args()

if sum(x is not None for x in [options.rotation,options.eulers,options.matrix,options.quaternion]) != 1:
  parser.error('more than one rotation specified.')

if options.quaternion is not None:
  eulers = damask.Rotation.fromQuaternion(np.array(options.quaternion)).asEulers(degrees=True)
if options.rotation is not None:
  eulers = damask.Rotation.fromAxisAngle(np.array(options.rotation,degrees=options.degrees)).asEulers(degrees=True)
if options.matrix is not None:
  eulers = damask.Rotation.fromMatrix(np.array(options.Matrix)).asEulers(degrees=True)
if options.eulers is not None:
  eulers = damask.Rotation.fromEulers(np.array(options.eulers),degrees=options.degrees).asEulers(degrees=True)

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
  
  spacing = geom.get_size()/geom.get_grid()

  newGrainID = options.fill if options.fill is not None else np.nanmax(microstructure)+1
  microstructure = ndimage.rotate(microstructure,eulers[2],(0,1),order=0,
                                  prefilter=False,output=microstructure.dtype,cval=newGrainID) # rotation around Z
  microstructure = ndimage.rotate(microstructure,eulers[1],(1,2),order=0,
                                  prefilter=False,output=microstructure.dtype,cval=newGrainID) # rotation around X
  microstructure = ndimage.rotate(microstructure,eulers[0],(0,1),order=0,
                                  prefilter=False,output=microstructure.dtype,cval=newGrainID) # rotation around Z
  
  geom.microstructure = microstructure
  geom.set_size(microstructure.shape*spacing)
  geom.add_comment(scriptID + ' ' + ' '.join(sys.argv[1:]))
  
  damask.util.croak(geom)
  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
