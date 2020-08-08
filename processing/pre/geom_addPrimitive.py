#!/usr/bin/env python3

import os
from optparse import OptionParser

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
                  help='i,j,k exponents for axes: '+
                  '0 gives octahedron (|x|^(2^0) + |y|^(2^0) + |z|^(2^0) < 1), '+
                  '1 gives a sphere (|x|^(2^1) + |y|^(2^1) + |z|^(2^1) < 1), '+
                  'large values produce boxes, negative turn concave.')
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

print('This script is broken and does nothing')
