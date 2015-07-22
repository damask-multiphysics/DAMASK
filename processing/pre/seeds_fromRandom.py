#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math,random
import numpy as np
import damask
from optparse import OptionParser,OptionGroup
from scipy import spatial


scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

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
parser.add_option('-w', '--weights', dest='weights', action='store_true',
                  help = 'assign random weigts (Gaussian Distribution) to seed points for laguerre tessellation [%default]')
parser.add_option('-m', '--microstructure', dest='microstructure', type='int',
                  help='first microstructure index [%default]', metavar='int')
parser.add_option('-s','--selective', dest='selective', action='store_true',
                  help = 'selective picking of seed points from random seed points [%default]')

group = OptionGroup(parser, "Laguerre Tessellation Options",
                   "Parameters determining shape of weight distribution of seed points "
                   )
group.add_option('--mean', dest='mean', type='float', metavar='float', \
                  help='mean of Gaussian Distribution for weights [%default]')
group.add_option('--sigma', dest='sigma', type='float', metavar='float', \
                  help='standard deviation of Gaussian Distribution for weights [%default]')
parser.add_option_group(group)

group = OptionGroup(parser, "Selective Seeding Options",
                    "More uniform distribution of seed points using Mitchell\'s Best Candidate Algorithm"
                   )
group.add_option('--distance', dest='bestDistance', type='float', metavar='float', \
                  help='minimum distance to the next neighbor  [%default]')
group.add_option('--numCandidates', dest='numCandidates', type='int', metavar='int', \
                  help='maximum number of point to consider for initial random points generation  [%default]')    
parser.add_option_group(group)

parser.set_defaults(randomSeed = None)
parser.set_defaults(grid = (16,16,16))
parser.set_defaults(N = 20)
parser.set_defaults(weights=False)
parser.set_defaults(mean = 0.0)
parser.set_defaults(sigma = 1.0)
parser.set_defaults(microstructure = 1)
parser.set_defaults(selective = False)
parser.set_defaults(bestDistance = 0.2)
parser.set_defaults(numCandidates = 10)



(options,filename) = parser.parse_args()
options.grid = np.array(options.grid)

labels = "1_coords\t2_coords\t3_coords\tphi1\tPhi\tphi2\tmicrostructure"

# ------------------------------------------ Functions Definitions ---------------------------------

def kdtree_search(xyz, point) :
    dist, index = spatial.cKDTree(xyz).query(np.array(point))
    return dist
    
def generatePoint() :
    return np.array([random.uniform(0,float(options.grid[0])/float(max(options.grid))), \
                     random.uniform(0,float(options.grid[1])/float(max(options.grid))), \
                     random.uniform(0,float(options.grid[2])/float(max(options.grid)))])


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
randomSeed = int(os.urandom(4).encode('hex'), 16)  if options.randomSeed == None else options.randomSeed
np.random.seed(randomSeed)                                                             # init random generators
random.seed(randomSeed)

grainEuler = np.random.rand(3,options.N)                                                # create random Euler triplets
grainEuler[0,:] *= 360.0                                                                # phi_1    is uniformly distributed
grainEuler[1,:] = np.arccos(2*grainEuler[1,:]-1)*180.0/math.pi                          # cos(Phi) is uniformly distributed
grainEuler[2,:] *= 360.0                                                                # phi_2    is uniformly distributed

microstructure=np.arange(options.microstructure,options.microstructure+options.N).reshape(1,options.N)

if options.selective == False :
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
  table = np.transpose(np.concatenate((seeds,grainEuler,microstructure),axis = 0))
else :
  samples = generatePoint().reshape(1,3)

  while  samples.shape[0] < options.N :
     bestDistance  = options.bestDistance
     for i in xrange(options.numCandidates) :
       c = generatePoint()
       d = kdtree_search(samples, c)
       if (d > bestDistance) :
        bestDistance = d
        bestCandidate = c
     if kdtree_search(samples,bestCandidate) != 0.0 :
        samples = np.append(samples,bestCandidate.reshape(1,3),axis=0)
     else :
        continue
  table = np.transpose(np.concatenate((samples.T,grainEuler,microstructure),axis = 0))

if options.weights :
   weight = np.random.normal(loc=options.mean, scale=options.sigma, size=options.N)
   table = np.append(table, weight.reshape(options.N,1), axis=1)
   labels += "\tweight"

# -------------------------------------- Write Data --------------------------------------------------

header = ["5\theader",
          scriptID + " " + " ".join(sys.argv[1:]),
          "grid\ta {}\tb {}\tc {}".format(options.grid[0],options.grid[1],options.grid[2]),
          "microstructures\t{}".format(options.N),
          "randomSeed\t{}".format(randomSeed),
          "%s"%labels,
         ]

for line in header:
  file['output'].write(line+"\n")
np.savetxt(file['output'], table, fmt='%10.6f', delimiter='\t')
