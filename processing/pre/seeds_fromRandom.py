#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math,random
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = scriptID.split()[1][:-3]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [options]', description = """
Distribute given number of points randomly within the three-dimensional cube [0.0,0.0,0.0]--[1.0,1.0,1.0].
Reports positions with random crystal orientations in seeds file format to STDOUT.

""", version = scriptID)

parser.add_option('-N', dest='N', type='int', metavar='int', \
                  help='number of seed points to distribute [%default]')
parser.add_option('-g','--grid', dest='grid', type='int', nargs=3, metavar='int int int', \
                  help='min a,b,c grid of hexahedral box %default')
parser.add_option('-r', '--rnd', dest='randomSeed', type='int', metavar='int', \
                  help='seed of random number generator [%default]')

parser.set_defaults(randomSeed = None)
parser.set_defaults(grid = (16,16,16))
parser.set_defaults(N = 20)

(options, extras) = parser.parse_args()

Npoints = reduce(lambda x, y: x * y, options.grid)
if 0 in options.grid: 
  file['croak'].write('invalid grid a b c.\n')
  sys.exit()
if options.N > Npoints: 
  sys.stderr.write('Warning: more seeds than grid points at minimum resolution.\n')
  options.N = Npoints
if options.randomSeed == None:
  options.randomSeed = int(os.urandom(4).encode('hex'), 16)

seeds = np.zeros((3,options.N),float)
np.random.seed(options.randomSeed)

grainEuler = np.random.rand(3,options.N)
grainEuler[0,:] *= 360.0
grainEuler[1,:] = np.arccos(2*grainEuler[1,:]-1)*180.0/math.pi
grainEuler[2,:] *= 360.0

seedpoint = np.random.permutation(Npoints)[:options.N]
seeds[0,:]=(np.mod(seedpoint                                   ,options.grid[0])\
                                                            +np.random.random())/options.grid[0]
seeds[1,:]=(np.mod(seedpoint//                 options.grid[0] ,options.grid[1])\
                                                            +np.random.random())/options.grid[1]
seeds[2,:]=(np.mod(seedpoint//(options.grid[1]*options.grid[0]),options.grid[2])\
                                                            +np.random.random())/options.grid[2]

sys.stdout.write("5\theader\n")
sys.stdout.write(scriptID + " " + " ".join(sys.argv[1:])+"\n")
sys.stdout.write("grid\ta {}\tb {}\tc {}\n".format(options.grid[0],options.grid[1],options.grid[2]))
sys.stdout.write("microstructures\t{}\n".format(options.N))
sys.stdout.write("randomSeed\t{}\n".format(options.randomSeed))
sys.stdout.write("x\ty\tz\tphi1\tPhi\tphi2\n")

np.savetxt(sys.stdout,np.transpose(np.concatenate((seeds,grainEuler),axis = 0)),fmt='%10.6f',delimiter='\t')
