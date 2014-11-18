#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math,random
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

(options,filename) = parser.parse_args()
options.grid = np.array(options.grid)

# ------------------------------------------ setup file handle -------------------------------------
if filename == []:
  file = {'output':sys.stdout, 'croak':sys.stderr}
else:
  file = {'output':open(filename[0],'w'), 'croak':sys.stderr}

gridSize = options.grid.prod()
if gridSize == 0:
  file['croak'].write('zero grid dimension for %s.\n'%(', '.join([['a','b','c'][x] for x in np.where(options.grid == 0)[0]])))
  sys.exit()
if options.N > gridSize: 
  file['croak'].write('accommodating only %i seeds on grid.\n'%gridSize)
  options.N = gridSize
if options.randomSeed == None:
  options.randomSeed = int(os.urandom(4).encode('hex'), 16)
np.random.seed(options.randomSeed)                                                      # init random generators
random.seed(options.randomSeed)

grainEuler = np.random.rand(3,options.N)                                                # create random Euler triplets
grainEuler[0,:] *= 360.0                                                                # phi_1    is uniformly distributed
grainEuler[1,:] = np.arccos(2*grainEuler[1,:]-1)*180.0/math.pi                          # cos(Phi) is uniformly distributed
grainEuler[2,:] *= 360.0                                                                # phi_2    is uniformly distributed

seedpoints = -np.ones(options.N,dtype='int')                                            # init grid positions of seed points

if options.N * 1024 < gridSize:                                                         # heuristic limit for random search
  i = 0
  while i < options.N:                                                                  # until all (unique) points determined
    p = np.random.randint(gridSize)                                                     # pick a location
    if p not in seedpoints:                                                             # not yet taken?
      seedpoints[i] = p                                                                 # take it
      i += 1                                                                            # advance stepper
else:
  seedpoints = np.array(random.sample(range(gridSize),options.N))                       # create random permutation of all grid positions and choose first N

seeds = np.zeros((3,options.N),float)                                                   # init seed positions
seeds[0,:] = (np.mod(seedpoints                                   ,options.grid[0])\
             +np.random.random())/options.grid[0]
seeds[1,:] = (np.mod(seedpoints//                 options.grid[0] ,options.grid[1])\
             +np.random.random())/options.grid[1]
seeds[2,:] = (np.mod(seedpoints//(options.grid[1]*options.grid[0]),options.grid[2])\
             +np.random.random())/options.grid[2]

header = ["5\theader",
          scriptID + " " + " ".join(sys.argv[1:]),
          "grid\ta {}\tb {}\tc {}".format(options.grid[0],options.grid[1],options.grid[2]),
          "microstructures\t{}".format(options.N),
          "randomSeed\t{}".format(options.randomSeed),
          "x\ty\tz\tphi1\tPhi\tphi2",
         ]

for line in header:
  file['output'].write(line+"\n")
np.savetxt(file['output'],np.transpose(np.concatenate((seeds,grainEuler),axis = 0)),fmt='%10.6f',delimiter='\t')
