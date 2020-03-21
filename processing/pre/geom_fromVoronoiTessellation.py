#!/usr/bin/env python3

import os
import sys
import multiprocessing
from io import StringIO
from optparse import OptionParser,OptionGroup

import numpy as np
from scipy import spatial

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


def Laguerre_tessellation(grid, seeds, grains, size, periodic, weights, cpus):

  def findClosestSeed(fargs):
    point, seeds, myWeights = fargs
    tmp = np.repeat(point.reshape(3,1), len(seeds), axis=1).T
    dist = np.sum((tmp - seeds)**2,axis=1) -myWeights
    return np.argmin(dist)                                                                          # seed point closest to point

  if periodic:
    weights_p = np.tile(weights,27).flatten(order='F')                                              # Laguerre weights (1,2,3,1,2,3,...,1,2,3)
    seeds_p = np.vstack((seeds  +np.array([size[0],0.,0.]),seeds,  seeds  +np.array([size[0],0.,0.])))
    seeds_p = np.vstack((seeds_p+np.array([0.,size[1],0.]),seeds_p,seeds_p+np.array([0.,size[1],0.])))
    seeds_p = np.vstack((seeds_p+np.array([0.,0.,size[2]]),seeds_p,seeds_p+np.array([0.,0.,size[2]])))
  else:
    weights_p = weights.flatten()
    seeds_p   = seeds

  arguments = [[arg,seeds_p,weights_p] for arg in list(grid)]

  if cpus > 1:                                                                                      # use multithreading
    pool = multiprocessing.Pool(processes = cpus)                                                   # initialize workers
    result = pool.map_async(findClosestSeed, arguments)                                             # evaluate function in parallel
    pool.close()
    pool.join()
    closestSeeds = np.array(result.get()).flatten()
  else:
    closestSeeds = np.zeros(len(arguments),dtype='i')
    for i,arg in enumerate(arguments):
      closestSeeds[i] = findClosestSeed(arg)

  return grains[closestSeeds%seeds.shape[0]]


def Voronoi_tessellation(grid, seeds, grains, size, periodic = True):

  KDTree = spatial.cKDTree(seeds,boxsize=size) if periodic else spatial.cKDTree(seeds)
  devNull,closestSeeds = KDTree.query(grid)
  return grains[closestSeeds]


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog option(s) [seedfile(s)]', description = """
Generate geometry description and material configuration by tessellation of given seeds file.

""", version = scriptID)


group = OptionGroup(parser, "Tessellation","")

group.add_option('-l',
                 '--laguerre',
                 dest = 'laguerre',
                 action = 'store_true',
                 help = 'use Laguerre (weighted Voronoi) tessellation')
group.add_option('--cpus',
                 dest = 'cpus',
                 type = 'int', metavar = 'int',
                 help = 'number of parallel processes to use for Laguerre tessellation [%default]')
group.add_option('--nonperiodic',
                 dest = 'periodic',
                 action = 'store_false',
                 help = 'nonperiodic tessellation')

parser.add_option_group(group)

group = OptionGroup(parser, "Geometry","")

group.add_option('-g',
                 '--grid',
                 dest = 'grid',
                 type = 'int', nargs = 3, metavar = ' '.join(['int']*3),
                 help = 'a,b,c grid of hexahedral box')
group.add_option('-s',
                 '--size',
                 dest = 'size',
                 type = 'float', nargs = 3, metavar=' '.join(['float']*3),
                 help = 'x,y,z size of hexahedral box [1.0 1.0 1.0]')
group.add_option('-o',
                 '--origin',
                 dest = 'origin',
                 type = 'float', nargs = 3, metavar=' '.join(['float']*3),
                 help = 'origin of grid [0.0 0.0 0.0]')

parser.add_option_group(group)

group = OptionGroup(parser, "Seeds","")

group.add_option('-p',
                 '--pos', '--seedposition',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'label of coordinates [%default]')
group.add_option('-w',
                 '--weight',
                 dest = 'weight',
                 type = 'string', metavar = 'string',
                 help = 'label of weights [%default]')
group.add_option('-m',
                 '--microstructure',
                 dest = 'microstructure',
                 type = 'string', metavar = 'string',
                 help = 'label of microstructures [%default]')
group.add_option('-e',
                 '--eulers',
                 dest = 'eulers',
                 type = 'string', metavar = 'string',
                 help = 'label of Euler angles [%default]')
group.add_option('--axes',
                 dest = 'axes',
                 type = 'string', nargs = 3, metavar = ' '.join(['string']*3),
                 help = 'orientation coordinate frame in terms of position coordinate frame')

parser.add_option_group(group)

group = OptionGroup(parser, "Configuration","")

group.add_option('--without-config',
                 dest = 'config',
                 action = 'store_false',
                 help = 'omit material configuration header')
group.add_option('--homogenization',
                 dest = 'homogenization',
                 type = 'int', metavar = 'int',
                 help = 'homogenization index to be used [%default]')
group.add_option('--phase',
                 dest = 'phase',
                 type = 'int', metavar = 'int',
                 help = 'phase index to be used [%default]')

parser.add_option_group(group)

parser.set_defaults(pos            = 'pos',
                    weight         = 'weight',
                    microstructure = 'microstructure',
                    eulers         = 'euler',
                    homogenization = 1,
                    phase          = 1,
                    cpus           = 2,
                    laguerre       = False,
                    periodic       = True,
                    config         = True,
                  )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)

    size   = np.ones(3)
    origin = np.zeros(3)
    for line in table.comments:
        items = line.lower().strip().split()
        key = items[0] if items else ''
        if   key == 'grid':
            grid   = np.array([  int(dict(zip(items[1::2],items[2::2]))[i]) for i in ['a','b','c']])
        elif key == 'size':
            size   = np.array([float(dict(zip(items[1::2],items[2::2]))[i]) for i in ['x','y','z']])
        elif key == 'origin':
            origin = np.array([float(dict(zip(items[1::2],items[2::2]))[i]) for i in ['x','y','z']])
    if options.grid:   grid   = np.array(options.grid)
    if options.size:   size   = np.array(options.size)
    if options.origin: origin = np.array(options.origin)


    seeds  = table.get(options.pos)

    grains = table.get(options.microstructure) if options.microstructure in table.labels else np.arange(len(seeds))+1
    grainIDs  = np.unique(grains).astype('i')
    NgrainIDs = len(grainIDs)

    if options.eulers in table.labels:
        eulers = table.get(options.eulers)

    coords = damask.grid_filters.cell_coord0(grid,size,-origin).reshape(-1,3,order='F')

    if options.laguerre:
        indices = Laguerre_tessellation(coords,seeds,grains,size,options.periodic,
                                        table.get(options.weight),options.cpus)
    else:
        indices = Voronoi_tessellation (coords,seeds,grains,size,options.periodic)

    config_header = []
    if options.config:

        if options.eulers in table.labels:
            config_header += ['<texture>']
            for ID in grainIDs:
                eulerID = np.nonzero(grains == ID)[0][0]                                            # find first occurrence of this grain id
                config_header += ['[Grain{}]'.format(ID),
                                  '(gauss)\tphi1 {:.2f}\tPhi {:.2f}\tphi2 {:.2f}'.format(*eulers[eulerID])
                                 ]
                if options.axes: config_header += ['axes\t{} {} {}'.format(*options.axes)]

        config_header += ['<microstructure>']
        for ID in grainIDs:
            config_header += ['[Grain{}]'.format(ID),
                              '(constituent)\tphase {}\ttexture {}\tfraction 1.0'.format(options.phase,ID)
                             ]

        config_header += ['<!skip>']

    header = [scriptID + ' ' + ' '.join(sys.argv[1:])]\
           + config_header
    geom = damask.Geom(indices.reshape(grid),size,origin,
                       homogenization=options.homogenization,comments=header)
    damask.util.croak(geom)

    geom.to_file(sys.stdout if name is None else os.path.splitext(name)[0]+'.geom',pack=False)
