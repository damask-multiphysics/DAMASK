#!/usr/bin/env python2
# -*- coding: UTF-8 no BOM -*-

import os,sys,math,random
import numpy as np
import damask
from optparse import OptionParser,OptionGroup
from scipy import spatial


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# ------------------------------------------ aux functions ---------------------------------

def kdtree_search(cloud, queryPoints):
  """find distances to nearest neighbor among cloud (N,d) for each of the queryPoints (n,d)"""
  n = queryPoints.shape[0]
  distances = np.zeros(n,dtype=float)
  tree = spatial.cKDTree(cloud)
  
  for i in xrange(n):
    distances[i], index = tree.query(queryPoints[i])

  return distances
    
# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [options]', description = """
Distribute given number of points randomly within the three-dimensional cube [0.0,0.0,0.0]--[1.0,1.0,1.0].
Reports positions with random crystal orientations in seeds file format to STDOUT.

""", version = scriptID)

parser.add_option('-N',
                  dest = 'N',
                  type = 'int', metavar = 'int',
                  help = 'number of seed points [%default]')
parser.add_option('-g',
                  '--grid',
                  dest = 'grid',
                  type = 'int', nargs = 3, metavar = 'int int int',
                  help='min a,b,c grid of hexahedral box %default')
parser.add_option('-m',
                  '--microstructure',
                  dest = 'microstructure',
                  type = 'int', metavar = 'int',
                  help = 'first microstructure index [%default]')
parser.add_option('-r',
                  '--rnd',
                  dest = 'randomSeed', type = 'int', metavar = 'int',
                  help = 'seed of random number generator [%default]')
parser.add_option('--format',
                  dest = 'format', type = 'string', metavar = 'string',
                  help = 'output number format [auto]')

group = OptionGroup(parser, "Laguerre Tessellation",
                   "Parameters determining shape of weight distribution of seed points"
                   )
group.add_option( '-w',
                  '--weights',
                  action = 'store_true',
                  dest   = 'weights',
                  help   = 'assign random weights to seed points for Laguerre tessellation [%default]')
group.add_option( '--max',
                  dest = 'max',
                  type = 'float', metavar = 'float',
                  help = 'max of uniform distribution for weights [%default]')
group.add_option( '--mean',
                  dest = 'mean',
                  type = 'float', metavar = 'float',
                  help = 'mean of normal distribution for weights [%default]')
group.add_option( '--sigma',
                  dest = 'sigma',
                  type = 'float', metavar = 'float',
                  help='standard deviation of normal distribution for weights [%default]')
parser.add_option_group(group)

group = OptionGroup(parser, "Selective Seeding",
                    "More uniform distribution of seed points using Mitchell's Best Candidate Algorithm"
                   )
group.add_option( '-s',
                  '--selective',
                  action = 'store_true',
                  dest   = 'selective',
                  help   = 'selective picking of seed points from random seed points [%default]')
group.add_option( '-f',
                  '--force',
                  action = 'store_true',
                  dest   = 'force',
                  help   = 'try selective picking despite large seed point number [%default]')
group.add_option( '--distance',
                  dest = 'distance',
                  type = 'float', metavar = 'float',
                  help = 'minimum distance to next neighbor [%default]')
group.add_option( '--numCandidates',
                  dest = 'numCandidates',
                  type = 'int', metavar = 'int',
                  help = 'size of point group to select best distance from [%default]')    
parser.add_option_group(group)

parser.set_defaults(randomSeed = None,
                    grid = (16,16,16),
                    N = 20,
                    weights = False,
                    max = 0.0,
                    mean = 0.2,
                    sigma = 0.05,
                    microstructure = 1,
                    selective = False,
                    force = False,
                    distance = 0.2,
                    numCandidates = 10,
                    format = None,
                   )

(options,filenames) = parser.parse_args()

options.grid = np.array(options.grid)
gridSize = options.grid.prod()

if options.randomSeed is None: options.randomSeed = int(os.urandom(4).encode('hex'), 16)
np.random.seed(options.randomSeed)                                                                  # init random generators
random.seed(options.randomSeed)


# --- loop over output files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(outname = name,
                                    buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# --- sanity checks -------------------------------------------------------------------------

  remarks = []
  errors  = []
  if gridSize == 0:
    errors.append('zero grid dimension for {}.'.format(', '.join([['a','b','c'][x] for x in np.where(options.grid == 0)[0]])))
  if options.N > gridSize/10.: errors.append('seed count exceeds 0.1 of grid points.')
  if options.selective and 4./3.*math.pi*(options.distance/2.)**3*options.N > 0.5:
    (remarks if options.force else errors).append('maximum recommended seed point count for given distance is {}.{}'.
                             format(int(3./8./math.pi/(options.distance/2.)**3),'..'*options.force))

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    sys.exit()

# --- do work ------------------------------------------------------------------------------------
 
  grainEuler = np.random.rand(3,options.N)                                                          # create random Euler triplets
  grainEuler[0,:] *= 360.0                                                                          # phi_1    is uniformly distributed
  grainEuler[1,:] = np.degrees(np.arccos(2*grainEuler[1,:]-1))                                      # cos(Phi) is uniformly distributed
  grainEuler[2,:] *= 360.0                                                                          # phi_2    is uniformly distributed

  if not options.selective:

    seeds = np.zeros((3,options.N),dtype='d')                                                       # seed positions array
    gridpoints = random.sample(range(gridSize),options.N)                                           # choose first N from random permutation of grid positions

    seeds[0,:] = (np.mod(gridpoints                                   ,options.grid[0])\
                 +np.random.random(options.N))                        /options.grid[0]
    seeds[1,:] = (np.mod(gridpoints//                 options.grid[0] ,options.grid[1])\
                 +np.random.random(options.N))                        /options.grid[1]
    seeds[2,:] = (np.mod(gridpoints//(options.grid[1]*options.grid[0]),options.grid[2])\
                 +np.random.random(options.N))                        /options.grid[2]

  else:

    seeds = np.zeros((options.N,3),dtype=float)                                                     # seed positions array
    seeds[0] = np.random.random(3)*options.grid/max(options.grid)
    i = 1                                                                                           # start out with one given point
    if i%(options.N/100.) < 1: damask.util.croak('.',False)

    while i < options.N:
      candidates = np.random.random(options.numCandidates*3).reshape(options.numCandidates,3)
      distances  = kdtree_search(seeds[:i],candidates)
      best = distances.argmax()
      if distances[best] > options.distance:                                                        # require minimum separation
        seeds[i] = candidates[best]                                                                 # maximum separation to existing point cloud
        i += 1
        if i%(options.N/100.) < 1: damask.util.croak('.',False)

    damask.util.croak('')
    seeds = seeds.T                                                                                 # prepare shape for stacking

  if options.weights:
    weights = [np.random.uniform(low = 0, high = options.max, size = options.N)] if options.max > 0.0 \
         else [np.random.normal(loc = options.mean, scale = options.sigma, size = options.N)]
  else:
    weights = []
  seeds = np.transpose(np.vstack(tuple([seeds,
                                        grainEuler,
                                        np.arange(options.microstructure,
                                                  options.microstructure + options.N),
                                       ] + weights
                                 )))

# ------------------------------------------ assemble header ---------------------------------------

  table.info_clear()
  table.info_append([
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {}\tb {}\tc {}".format(*options.grid),
    "microstructures\t{}".format(options.N),
    "randomSeed\t{}".format(options.randomSeed),
    ])
  table.labels_clear()
  table.labels_append( ['{dim}_{label}'.format(dim = 1+k,label = 'pos')   for k in xrange(3)] +
                       ['{dim}_{label}'.format(dim = 1+k,label = 'euler') for k in xrange(3)] + 
                       ['microstructure'] +
                      (['weight'] if options.weights else []))
  table.head_write()
  table.output_flush()
  
# --- write seeds information ------------------------------------------------------------

  table.data = seeds
  table.data_writeArray(fmt = options.format)
    
# --- output finalization --------------------------------------------------------------------------

  table.close()                                                                                     # close ASCII table
