#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math,string
import numpy as np
import multiprocessing
from optparse import OptionParser
from scipy    import spatial
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]


def meshgrid2(*arrs):
  '''
  code inspired by http://stackoverflow.com/questions/1827489/numpy-meshgrid-in-3d
  '''
  arrs = tuple(reversed(arrs))
  arrs = tuple(arrs)
  lens = np.array(map(len, arrs))
  dim = len(arrs)
  ans = []
  for i, arr in enumerate(arrs):
    slc = np.ones(dim,'i')
    slc[i] = lens[i]
    arr2 = np.asarray(arr).reshape(slc)
    for j, sz in enumerate(lens):
      if j != i:
        arr2 = arr2.repeat(sz, axis=j)
   
    ans.insert(0,arr2)
  return tuple(ans)

def findClosestSeed(fargs):
    point, seeds, weightssquared = fargs
    tmp = np.repeat(point.reshape(3,1), len(seeds), axis=1).T
    dist = np.sum((tmp - seeds)*(tmp - seeds),axis=1) - weightssquared
    return np.argmin(dist)                                                                          # seed point closest to point


def laguerreTessellation(undeformed, coords, weights, grains, nonperiodic = False, cpus = 2):

    copies = \
      np.array([
                [  0, 0, 0 ],
               ]).astype(float) if nonperiodic else \
     np.array([
                [ -1,-1,-1 ],
                [  0,-1,-1 ],
                [  1,-1,-1 ],
                [ -1, 0,-1 ],
                [  0, 0,-1 ],
                [  1, 0,-1 ],
                [ -1, 1,-1 ],
                [  0, 1,-1 ],
                [  1, 1,-1 ],
                [ -1,-1, 0 ],
                [  0,-1, 0 ],
                [  1,-1, 0 ],
                [ -1, 0, 0 ],
                [  0, 0, 0 ],
                [  1, 0, 0 ],
                [ -1, 1, 0 ],
                [  0, 1, 0 ],
                [  1, 1, 0 ],
                [ -1,-1, 1 ],
                [  0,-1, 1 ],
                [  1,-1, 1 ],
                [ -1, 0, 1 ],
                [  0, 0, 1 ],
                [  1, 0, 1 ],
                [ -1, 1, 1 ],
                [  0, 1, 1 ],
                [  1, 1, 1 ],
               ]).astype(float)*info['size']

    squaredweights = np.power(np.tile(weights,len(copies)),2)                                       # Laguerre weights (squared, size N*n)

    for i,vec in enumerate(copies):                                                                 # periodic copies of seed points (size N*n)
      try: seeds = np.append(seeds, coords+vec, axis=0)
      except NameError: seeds = coords+vec   

    if all(squaredweights == 0.0):                                                                  # standard Voronoi (no weights, KD tree)
      myKDTree = spatial.cKDTree(seeds)
      devNull,closestSeeds = myKDTree.query(undeformed)
    else:
      damask.util.croak('...using {} cpu{}'.format(options.cpus, 's' if options.cpus > 1 else ''))
      arguments = [[arg] + [seeds,squaredweights] for arg in list(undeformed)]

      if cpus > 1:                                                                                  # use multithreading
        pool = multiprocessing.Pool(processes = cpus)                                               # initialize workers
        result = pool.map_async(findClosestSeed, arguments)                                         # evaluate function in parallel
        pool.close()
        pool.join()
        closestSeeds = np.array(result.get()).flatten()
      else:
        closestSeeds = np.zeros(len(arguments),dtype='i')
        for i,arg in enumerate(arguments):
          closestSeeds[i] = findClosestSeed(arg)

    return grains[closestSeeds%coords.shape[0]]                                                   # closestSeed is modulo number of original seed points (i.e. excluding periodic copies)

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Generate geometry description and material configuration by standard Voronoi tessellation of given seeds file.

""", version = scriptID)

parser.add_option('-g', '--grid',
                  dest = 'grid',
                  type = 'int', nargs = 3, metavar = ' '.join(['int']*3),
                  help = 'a,b,c grid of hexahedral box [from seeds file]')
parser.add_option('-s', '--size',
                  dest = 'size',
                  type = 'float', nargs = 3, metavar=' '.join(['float']*3),
                  help = 'x,y,z size of hexahedral box [from seeds file or 1.0 along largest grid point number]')
parser.add_option('-o', '--origin',
                  dest = 'origin',
                  type = 'float', nargs = 3, metavar=' '.join(['float']*3),
                  help = 'offset from old to new origin of grid')
parser.add_option('-p', '--position',
                  dest = 'position',
                  type = 'string', metavar = 'string',
                  help = 'column label for seed positions [%default]')
parser.add_option('-w', '--weight',
                  dest = 'weight',
                  type = 'string', metavar = 'string',
                  help = 'column label for seed weights [%default]')
parser.add_option('-m', '--microstructure',
                  dest = 'microstructure',
                  type = 'string', metavar = 'string',
                  help = 'column label for seed microstructures [%default]')
parser.add_option('-e', '--eulers',
                  dest = 'eulers',
                  type = 'string', metavar = 'string',
                  help = 'column label for seed Euler angles [%default]')
parser.add_option('--axes',
                  dest = 'axes',
                  type = 'string', nargs = 3, metavar = ' '.join(['string']*3),
                  help = 'orientation coordinate frame in terms of position coordinate frame [same]')
parser.add_option('--homogenization',
                  dest = 'homogenization',
                  type = 'int', metavar = 'int',
                  help = 'homogenization index to be used [%default]')
parser.add_option('--crystallite',
                  dest = 'crystallite',
                  type = 'int', metavar = 'int',
                  help = 'crystallite index to be used [%default]')
parser.add_option('--phase',
                  dest = 'phase',
                  type = 'int', metavar = 'int',
                  help = 'phase index to be used [%default]')
parser.add_option('-r', '--rnd',
                  dest = 'randomSeed',
                  type = 'int', metavar='int',
                  help = 'seed of random number generator for second phase distribution [%default]')
parser.add_option('--secondphase',
                  dest = 'secondphase',
                  type = 'float', metavar= 'float',
                  help = 'volume fraction of randomly distribute second phase [%default]')
parser.add_option('-l', '--laguerre',
                  dest = 'laguerre',
                  action = 'store_true',
                  help = 'use Laguerre (weighted Voronoi) tessellation [%default]')
parser.add_option('--cpus',
                  dest = 'cpus',
                  type = 'int', metavar = 'int',
                  help = 'number of parallel processes to use for Laguerre tessellation [%default]')
parser.add_option('--nonperiodic',
                  dest = 'nonperiodic',
                  action = 'store_true',
                  help = 'use nonperiodic tessellation [%default]')

parser.set_defaults(grid   = None,
                    size   = None,
                    origin = None,
                    position       = 'pos',
                    weight         = 'weight',
                    microstructure = 'microstructure',
                    eulers         = 'eulerangles',
                    homogenization = 1,
                    crystallite    = 1,
                    phase          = 1,
                    secondphase    = 0.0,
                    cpus           = 2,
                    laguerre       = False,
                    nonperiodic    = False,
                    randomSeed     = None,
                  )
(options,filenames) = parser.parse_args()

if options.secondphase > 1.0 or options.secondphase < 0.0:
 parser.error('volume fraction of second phase ({}) out of bounds.'.format(options.secondphase))

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              outname = os.path.splitext(name)[-2]+'.geom' if name else name,
                              buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# --- read header ----------------------------------------------------------------------------

  table.head_read()
  info,extra_header = table.head_getGeom()
  
  if options.grid   != None: info['grid']   = options.grid
  if options.size   != None: info['size']   = options.size
  if options.origin != None: info['origin'] = options.origin
  
# ------------------------------------------ sanity checks ---------------------------------------  

  remarks = []
  errors = []
  labels = []
  
  hasGrains  = table.label_dimension(options.microstructure) == 1
  hasEulers  = table.label_dimension(options.eulers) == 3
  hasWeights = table.label_dimension(options.weight) == 1 and options.laguerre

  if     np.any(info['grid'] < 1):   errors.append('invalid grid a b c.')
  if     np.any(info['size'] <= 0.0) \
     and np.all(info['grid'] < 1):   errors.append('invalid size x y z.')
  else:
    for i in xrange(3):
      if info['size'][i] <= 0.0:                                                                      # any invalid size?
        info['size'][i] = float(info['grid'][i])/max(info['grid'])                                    # normalize to grid
        remarks.append('rescaling size {} to {}...'.format({0:'x',1:'y',2:'z'}[i],info['size'][i]))

  if table.label_dimension(options.position) != 3:
    errors.append('position columns "{}" have dimension {}.'.format(options.position,
                                                                    table.label_dimension(options.position)))
  else:
    labels += [options.position]

  if not hasEulers:                        remarks.append('missing seed orientations...')
  else: labels += [options.eulers]
  if not hasGrains:                        remarks.append('missing seed microstructure indices...')
  else: labels += [options.microstructure]
  if options.laguerre and not hasWeights:  remarks.append('missing seed weights...')
  else: labels += [options.weight]

  if remarks != []: damask.util.croak(remarks)
  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss=True)
    continue

# ------------------------------------------ read seeds ---------------------------------------  
      
  table.data_readArray(labels)
  coords = table.data[:,table.label_index(options.position):table.label_index(options.position)+3]\
           * info['size']
  eulers = table.data[:,table.label_index(options.eulers  ):table.label_index(options.eulers  )+3]\
           if hasEulers  else np.zeros(3*len(coords))
  grains = table.data[:,table.label_index(options.microstructure)].astype('i')\
           if hasGrains  else 1+np.arange(len(coords))
  weights = table.data[:,table.label_index(options.weight)]\
           if hasWeights else np.zeros(len(coords))
  grainIDs = np.unique(grains).astype('i')
  NgrainIDs = len(grainIDs)

# --- tessellate microstructure ------------------------------------------------------------

  x = (np.arange(info['grid'][0])+0.5)*info['size'][0]/info['grid'][0]
  y = (np.arange(info['grid'][1])+0.5)*info['size'][1]/info['grid'][1]
  z = (np.arange(info['grid'][2])+0.5)*info['size'][2]/info['grid'][2]
  
  damask.util.croak('tessellating...')

  grid = np.vstack(meshgrid2(x, y, z)).reshape(3,-1).T
  indices = laguerreTessellation(grid, coords, weights, grains, options.nonperiodic, options.cpus)
    
# --- write header ---------------------------------------------------------------------------------

  grainIDs = np.intersect1d(grainIDs,indices)
  info['microstructures'] = len(grainIDs)

  if info['homogenization'] == 0: info['homogenization'] = options.homogenization
  
  damask.util.croak(['grid     a b c:  %s'%(' x '.join(map(str,info['grid']))),
               'size     x y z:  %s'%(' x '.join(map(str,info['size']))),
               'origin   x y z:  %s'%(' : '.join(map(str,info['origin']))),
               'homogenization:  %i'%info['homogenization'],
               'microstructures: %i%s'%(info['microstructures'],
                                        (' out of %i'%NgrainIDs if NgrainIDs != info['microstructures'] else '')),
              ])

  config_header = []
  formatwidth = 1+int(math.log10(info['microstructures']))

  phase = options.phase * np.ones(info['microstructures'],'i')
  if int(options.secondphase*info['microstructures']) > 0:
    phase[0:int(options.secondphase*info['microstructures'])] += 1
    randomSeed = int(os.urandom(4).encode('hex'), 16)  if options.randomSeed == None \
                                                       else options.randomSeed                      # random seed for second phase
    np.random.seed(randomSeed)
    np.random.shuffle(phase)
    config_header += ['# random seed (phase shuffling): {}'.format(randomSeed)]

  config_header += ['<microstructure>']
  for i,ID in enumerate(grainIDs):
    config_header += ['[Grain%s]'%(str(ID).zfill(formatwidth)),
                      'crystallite %i'%options.crystallite,
                      '(constituent)\tphase %i\ttexture %s\tfraction 1.0'%(phase[i],str(ID).rjust(formatwidth)),
                     ]
  if hasEulers:
    config_header += ['<texture>']
    for ID in grainIDs:
      eulerID = np.nonzero(grains == ID)[0][0]                                                      # find first occurrence of this grain id
      config_header += ['[Grain%s]'%(str(ID).zfill(formatwidth)),
                        '(gauss)\tphi1 %g\tPhi %g\tphi2 %g\tscatter 0.0\tfraction 1.0'%tuple(eulers[eulerID])
                       ]
      if options.axes != None: config_header.append('axes\t%s %s %s'%tuple(options.axes))
  
  table.labels_clear()
  table.info_clear()
  table.info_append([
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=info['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=info['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=info['origin']),
    "homogenization\t{homog}".format(homog=info['homogenization']),
    "microstructures\t{microstructures}".format(microstructures=info['microstructures']),
    config_header,
    ])
  table.head_write()
      
# --- write microstructure information ------------------------------------------------------------

  table.data = indices.reshape(info['grid'][1]*info['grid'][2],info['grid'][0])
  table.data_writeArray('%%%ii'%(formatwidth),delimiter=' ')
    
#--- output finalization --------------------------------------------------------------------------

  table.close()  
