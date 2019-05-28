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
These objects can be boxes, cylinders, or ellipsoids.

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


if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)
  
  geom = damask.Geom.from_file(StringIO(''.join(sys.stdin.read())) if name is None else name)
  grid = geom.get_grid()
  size = geom.get_size()
  microstructure = geom.get_microstructure()

  # scale to box of size [1.0,1.0,1.0]
  if options.realspace:
    center = (np.array(options.center) - geom.get_origin())/size
    r = np.array(options.dimension)/size/2.0
  else:
    center = (np.array(options.center) + 0.5)/grid
    r = np.array(options.dimension)/grid/2.0
    
  if np.any(center<0.0) or np.any(center>=1.0): print('error')

  offset = np.ones(3)*0.5 if options.periodic else center
  mask = np.full(grid,False)
  # High exponents can cause underflow & overflow - okay here, just compare to 1, so +infinity and 0 are fine
  np.seterr(over='ignore', under='ignore')

  e = np.array(options.exponent)
  for x in range(grid[0]):
    for y in range(grid[1]):
      for z in range(grid[2]):
        coords = np.array([x+0.5,y+0.5,z+0.5])/grid
        mask[x,y,z] = np.sum(np.abs((rotation*(coords-offset))/r)**e) < 1
        
  if options.periodic:
    shift = ((offset-center)*grid).astype(int)
    mask = np.roll(mask,shift,(0,1,2))
    
    
  if options.inside: mask = np.logical_not(mask)
  fill = np.nanmax(microstructure)+1 if options.fill is None else options.fill
  microstructure = np.where(mask,microstructure,fill)

  damask.util.croak(geom.update(microstructure))
  geom.add_comments(scriptID + ' ' + ' '.join(sys.argv[1:]))
  
  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
