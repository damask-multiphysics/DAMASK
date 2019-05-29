#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

from scipy import ndimage
import numpy as np

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [geomfile(s)]', description = """
Rotates original microstructure and embeddeds it into buffer material.

""", version=scriptID)

parser.add_option('-r', '--rotation',
                  dest='rotation',
                  type = 'float', nargs = 4, metavar = ' '.join(['float']*4),
                  help = 'rotation given as axis and angle')
parser.add_option('-e', '--eulers',
                  dest = 'eulers',
                  type = 'float', nargs = 3, metavar = ' '.join(['float']*3),
                  help = 'rotation given as Euler angles')
parser.add_option('-d', '--degrees',
                  dest = 'degrees',
                  action = 'store_true',
                  help = 'Angles (Euler angles/axis angle) are given in degrees [%default]')
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
                  type = 'float', metavar = 'int',
                  help = 'background microstructure index, defaults to max microstructure index + 1')

parser.set_defaults(degrees = False)

(options, filenames) = parser.parse_args()

if [options.rotation,options.eulers,options.matrix,options.quaternion].count(None) < 3:
  parser.error('more than one rotation specified.')
if [options.rotation,options.eulers,options.matrix,options.quaternion].count(None) > 3:
  parser.error('no rotation specified.')

if options.quaternion is not None:
  rot = damask.Rotation.fromQuaternion(np.array(options.quaternion))                            # we might need P=+1 here, too...
if options.rotation is not None:
  rot = damask.Rotation.fromAxisAngle(np.array(options.rotation),degrees=options.degrees,P=+1)
if options.matrix is not None:
  rot = damask.Rotation.fromMatrix(np.array(options.Matrix))
if options.eulers is not None:
  rot = damask.Rotation.fromEulers(np.array(options.eulers),degrees=options.degrees)

eulers = rot.asEulers(degrees=True)

if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)

  geom = damask.Geom.from_file(StringIO(''.join(sys.stdin.read())) if name is None else name)
  microstructure = geom.get_microstructure()
  fill = np.nanmax(microstructure)+1 if options.fill is None else options.fill
  dtype = float if np.isnan(fill) or int(fill) != fill or microstructure.dtype==np.float else int

  # These rotations are always applied in the reference coordinate system, i.e. (z,x,z) not (z,x',z'')
  # this seems to be ok, see https://www.cs.utexas.edu/~theshark/courses/cs354/lectures/cs354-14.pdf
  microstructure = ndimage.rotate(microstructure,eulers[2],(0,1),order=0,
                                  prefilter=False,output=dtype,cval=fill)            # rotation around z
  microstructure = ndimage.rotate(microstructure,eulers[1],(1,2),order=0,
                                  prefilter=False,output=dtype,cval=fill)            # rotation around x
  microstructure = ndimage.rotate(microstructure,eulers[0],(0,1),order=0,
                                  prefilter=False,output=dtype,cval=fill)            # rotation around z
  
  damask.util.croak(geom.update(microstructure,rescale=True))
  geom.add_comments(scriptID + ' ' + ' '.join(sys.argv[1:]))
  
  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
