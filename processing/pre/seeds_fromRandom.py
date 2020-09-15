#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser,OptionGroup

import numpy as np
from scipy import spatial

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options', description = """
Distribute given number of points randomly within rectangular cuboid.
Reports positions with random crystal orientations in seeds file format to STDOUT.

""", version = scriptID)

parser.add_option('-N',
                  dest = 'N',
                  type = 'int', metavar = 'int',
                  help = 'number of seed points [%default]')
parser.add_option('-s',
                  '--size',
                  dest = 'size',
                  type = 'float', nargs = 3, metavar = 'float float float',
                  help='size x,y,z of unit cube to fill %default')
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
group.add_option( '--selective',
                  action = 'store_true',
                  dest   = 'selective',
                  help   = 'selective picking of seed points from random seed points')
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
                    size = (1.0,1.0,1.0),
                    N = 20,
                    weights = False,
                    max = 0.0,
                    mean = 0.2,
                    sigma = 0.05,
                    microstructure = 1,
                    selective = False,
                    distance = 0.2,
                    numCandidates = 10,
                   )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

size = np.array(options.size)
grid = np.array(options.grid)
np.random.seed(int(os.urandom(4).hex(),16) if options.randomSeed is None else options.randomSeed)

for name in filenames:
    damask.util.report(scriptName,name)

    if options.N > np.prod(grid):
        damask.util.croak('More seeds than grid positions.')
        sys.exit()
    if options.selective and options.distance < min(size/grid):
        damask.util.croak('Distance must be larger than grid spacing.')
        sys.exit()
    if options.selective and options.distance**3*options.N > 0.5*np.prod(size):
        damask.util.croak('Number of seeds for given size and distance should be < {}.'\
                          .format(int(0.5*np.prod(size)/options.distance**3)))

    eulers = np.random.rand(options.N,3)                                                            # create random Euler triplets
    eulers[:,0] *= 360.0                                                                            # phi_1    is uniformly distributed
    eulers[:,1] = np.degrees(np.arccos(2*eulers[:,1]-1.0))                                          # cos(Phi) is uniformly distributed
    eulers[:,2] *= 360.0                                                                            # phi_2    is uniformly distributed


    if not options.selective:
        coords = damask.grid_filters.cell_coord0(grid,size).reshape(-1,3,order='F')
        seeds = coords[np.random.choice(np.prod(grid), options.N, replace=False)] \
              + np.broadcast_to(size/grid,(options.N,3))*(np.random.rand(options.N,3)*.5-.25)       # wobble without leaving grid
    else:
        seeds = np.empty((options.N,3))
        seeds[0] = np.random.random(3) * size

        i = 1
        progress = damask.util._ProgressBar(options.N,'',50)
        while i < options.N:
            candidates = np.random.rand(options.numCandidates,3)*np.broadcast_to(size,(options.numCandidates,3))
            tree = spatial.cKDTree(seeds[:i])
            distances, dev_null = tree.query(candidates)
            best = distances.argmax()
            if distances[best] > options.distance:                                                  # require minimum separation
                seeds[i] = candidates[best]                                                         # maximum separation to existing point cloud
                i += 1
                progress.update(i)


    comments = [scriptID + ' ' + ' '.join(sys.argv[1:]),
                'grid\ta {}\tb {}\tc {}'.format(*grid),
                'size\tx {}\ty {}\tz {}'.format(*size),
                'randomSeed\t{}'.format(options.randomSeed),
               ]

    table = damask.Table(np.hstack((seeds,eulers)),{'pos':(3,),'euler':(3,)},comments)\
                  .add('microstructure',np.arange(options.microstructure,options.microstructure + options.N,dtype=int))

    if options.weights:
        weights = np.random.uniform(low = 0, high = options.max, size = options.N) if options.max > 0.0 \
             else np.random.normal(loc = options.mean, scale = options.sigma, size = options.N)
        table = table.add('weight',weights)

    table.save_ASCII(sys.stdout if name is None else name)
