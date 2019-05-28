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
Inserts a primitive geometric object at a given position.
Depending on the sign of the dimension parameters, these objects can be boxes, cylinders, or ellipsoids.

""", version = scriptID)

parser.add_option('-c', '--center',
                  dest='center', 
                  type='float', nargs = 3, metavar=' '.join(['float']*3),
                  help='a,b,c origin of primitive %default')
parser.add_option('-d', '--dimension',
                  dest='dimension',
                  type='float', nargs = 3, metavar=' '.join(['float']*3),
                  help='a,b,c extension of hexahedral box')
parser.add_option('-e', '--exponent',
                  dest='exponent',
                  type='float', nargs = 3, metavar=' '.join(['float']*3),
                  help='i,j,k exponents for axes - 0 gives octahedron (|x|^(2^0) + |y|^(2^0) + |z|^(2^0) < 1), \
                  1 gives a sphere (|x|^(2^1) + |y|^(2^1) + |z|^(2^1) < 1), \
                  large values produce boxes, negative turn concave.')
parser.add_option('-f', '--fill',
                  dest='fill',
                  type='float', metavar = 'int',
                  help='microstructure index to fill primitive, defaults to max microstructure index + 1')
parser.add_option('-q', '--quaternion',
                  dest='quaternion',
                  type='float', nargs = 4, metavar=' '.join(['float']*4),
                  help = 'rotation of primitive as quaternion')
parser.add_option('-a', '--angleaxis',
                  dest='angleaxis',
                  type=float, nargs = 4, metavar=' '.join(['float']*4),
                  help = 'axis and angle to rotate primitive')
parser.add_option(     '--degrees',
                  dest='degrees',
                  action='store_true',
                  help = 'angle is given in degrees')
parser.add_option(     '--nonperiodic',
                  dest='periodic',
                  action='store_false',
                  help = 'wrap around edges')
parser.add_option(     '--realspace',
                  dest='realspace',
                  action='store_true',
                  help = '-c and -d span [origin,origin+size] instead of [0,grid] coordinates')
parser.add_option(     '--invert',
                  dest='inside',
                  action='store_false',
                  help = 'invert the volume filled by the primitive (inside/outside)')

parser.set_defaults(center = (.0,.0,.0),
                    degrees = False,
                    exponent = (20,20,20), # box shape by default
                    periodic = True,
                    realspace = False,
                    inside = True,
                   )

(options, filenames) = parser.parse_args()

if options.dimension is None:
  parser.error('no dimension specified.') 
if [options.angleaxis,options.quaternion].count(None) == 0:
  parser.error('more than one rotation specified.')

if options.angleaxis is not None:
  rotation = damask.Rotation.fromAxisAngle(np.array(options.angleaxis),options.degrees,normalise=True)
elif options.quaternion is not None:
  rotation = damask.Rotation.fromQuaternion(options.quaternion)
else:
  rotation = damask.Rotation()

options.center    = np.array(options.center)
options.dimension = np.array(options.dimension)
# undo logarithmic sense of exponent and generate ellipsoids for negative dimensions (backward compatibility)
options.exponent  = np.where(np.array(options.dimension) > 0, np.power(2,options.exponent), 2)


if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)
  
  geom = damask.Geom.from_file(StringIO(''.join(sys.stdin.read())) if name is None else name)
  grid = geom.get_grid()
  size = geom.get_size()
  origin = geom.get_origin()
  microstructure = geom.get_microstructure()

  # coordinates given in real space, not (default) voxel space
  if options.realspace:
    options.center    -= origin
    options.center    *= grid / size
    options.dimension *= grid / size
  

  # change to coordinate space where the primitive is the unit sphere/cube/etc
  if options.periodic: # use padding to achieve periodicity
    (X, Y, Z) = np.meshgrid(np.arange(-grid[0]/2, (3*grid[0])/2, dtype=np.float32), # 50% padding on each side 
                            np.arange(-grid[1]/2, (3*grid[1])/2, dtype=np.float32), 
                            np.arange(-grid[2]/2, (3*grid[2])/2, dtype=np.float32),
                            indexing='ij')
    # Padding handling
    X = np.roll(np.roll(np.roll(X,
            -grid[0]//2, axis=0), 
            -grid[1]//2, axis=1), 
            -grid[2]//2, axis=2)
    Y = np.roll(np.roll(np.roll(Y,
            -grid[0]//2, axis=0), 
            -grid[1]//2, axis=1), 
            -grid[2]//2, axis=2)
    Z = np.roll(np.roll(np.roll(Z,
            -grid[0]//2, axis=0), 
            -grid[1]//2, axis=1), 
            -grid[2]//2, axis=2)
  else: # nonperiodic, much lighter on resources
    # change to coordinate space where the primitive is the unit sphere/cube/etc
    (X, Y, Z) = np.meshgrid(np.arange(0, grid[0], dtype=np.float32), 
                            np.arange(0, grid[1], dtype=np.float32), 
                            np.arange(0, grid[2], dtype=np.float32),
                            indexing='ij')
    
  # first by translating the center onto 0, 0.5 shifts the voxel origin onto the center of the voxel
  X -= options.center[0] - 0.5
  Y -= options.center[1] - 0.5
  Z -= options.center[2] - 0.5
  # and then by applying the rotation
  (X, Y, Z) = rotation * (X, Y, Z)
  # and finally by scaling (we don't worry about options.dimension being negative, np.abs occurs on the microstructure = np.where... line)
  X /= options.dimension[0] * 0.5
  Y /= options.dimension[1] * 0.5
  Z /= options.dimension[2] * 0.5
    
  fill = np.nanmax(microstructure)+1 if options.fill is None else options.fill

 # High exponents can cause underflow & overflow - loss of precision is okay here, we just compare it to 1, so +infinity and 0 are fine
  old_settings = np.seterr()
  np.seterr(over='ignore', under='ignore')
  
  if options.periodic: # use padding to achieve periodicity
    inside = np.zeros(grid, dtype=bool)
    for i in range(2):
      for j in range(2):
        for k in range(2):
          inside = inside | ( # Most of this is handling the padding
                np.abs(X[grid[0] * i : grid[0] * (i+1),
                         grid[1] * j : grid[1] * (j+1),
                         grid[2] * k : grid[2] * (k+1)])**options.exponent[0] +
                np.abs(Y[grid[0] * i : grid[0] * (i+1),
                         grid[1] * j : grid[1] * (j+1),
                         grid[2] * k : grid[2] * (k+1)])**options.exponent[1] +
                np.abs(Z[grid[0] * i : grid[0] * (i+1),
                         grid[1] * j : grid[1] * (j+1),
                         grid[2] * k : grid[2] * (k+1)])**options.exponent[2] <= 1.0)
    
    microstructure = np.where(inside,
                              fill           if options.inside else microstructure,
                              microstructure if options.inside else fill)   

  else: # nonperiodic, much lighter on resources
    microstructure = np.where(np.abs(X)**options.exponent[0] +
                              np.abs(Y)**options.exponent[1] +
                              np.abs(Z)**options.exponent[2] <= 1.0,
                              fill           if options.inside else microstructure,
                              microstructure if options.inside else fill)   

  damask.util.croak(geom.update(microstructure))
  geom.add_comments(scriptID + ' ' + ' '.join(sys.argv[1:]))
  
  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
