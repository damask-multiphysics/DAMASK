#!/usr/bin/env python3

import os
import sys 
import multiprocessing
from optparse import OptionParser,OptionGroup

import numpy as np
from scipy import spatial

import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


def meshgrid2(*arrs):
  """Code inspired by http://stackoverflow.com/questions/1827489/numpy-meshgrid-in-3d"""
  arrs = tuple(reversed(arrs))
  lens = np.array(list(map(len, arrs)))
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


def laguerreTessellation(undeformed, coords, weights, grains, periodic = True, cpus = 2):

  def findClosestSeed(fargs):
    point, seeds, myWeights = fargs
    tmp = np.repeat(point.reshape(3,1), len(seeds), axis=1).T
    dist = np.sum((tmp - seeds)**2,axis=1) -myWeights
    return np.argmin(dist)                                                                          # seed point closest to point
    
  copies = \
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
             ]).astype(float)*info['size'] if periodic else \
    np.array([
              [  0, 0, 0 ],
             ]).astype(float)

  repeatweights = np.tile(weights,len(copies)).flatten(order='F')                                   # Laguerre weights (1,2,3,1,2,3,...,1,2,3)
  for i,vec in enumerate(copies):                                                                   # periodic copies of seed points ...
    try: seeds = np.append(seeds, coords+vec, axis=0)                                               # ... (1+a,2+a,3+a,...,1+z,2+z,3+z)
    except NameError: seeds = coords+vec   

  if (repeatweights == 0.0).all():                                                                  # standard Voronoi (no weights, KD tree)
    myKDTree = spatial.cKDTree(seeds)
    devNull,closestSeeds = myKDTree.query(undeformed)
  else:
    damask.util.croak('...using {} cpu{}'.format(options.cpus, 's' if options.cpus > 1 else ''))
    arguments = [[arg,seeds,repeatweights] for arg in list(undeformed)]

    if cpus > 1:                                                                                    # use multithreading
      pool = multiprocessing.Pool(processes = cpus)                                                 # initialize workers
      result = pool.map_async(findClosestSeed, arguments)                                           # evaluate function in parallel
      pool.close()
      pool.join()
      closestSeeds = np.array(result.get()).flatten()
    else:
      closestSeeds = np.zeros(len(arguments),dtype='i')
      for i,arg in enumerate(arguments):
        closestSeeds[i] = findClosestSeed(arg)

# closestSeed is modulo number of original seed points (i.e. excluding periodic copies)
  return grains[closestSeeds%coords.shape[0]]                                                   


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
                 help = 'x,y,z size of hexahedral box')
group.add_option('-o',
                 '--origin',
                 dest = 'origin',
                 type = 'float', nargs = 3, metavar=' '.join(['float']*3),
                 help = 'origin of grid')
group.add_option('--nonnormalized',
                 dest = 'normalized',
                 action = 'store_false',
                 help = 'seed coordinates are not normalized to a unit cube')

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
                    normalized     = True,
                    config         = True,
                  )
(options,filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  table = damask.ASCIItable(name = name, readonly = True)
  damask.util.report(scriptName,name)

# --- read header ----------------------------------------------------------------------------

  table.head_read()
  info,extra_header = table.head_getGeom()
  
  if options.grid   is not None: info['grid']   = options.grid
  if options.size   is not None: info['size']   = options.size
  if options.origin is not None: info['origin'] = options.origin
  
# ------------------------------------------ sanity checks ---------------------------------------  

  remarks = []
  errors = []
  labels = []
  
  hasGrains  = table.label_dimension(options.microstructure) == 1
  hasEulers  = table.label_dimension(options.eulers) == 3
  hasWeights = table.label_dimension(options.weight) == 1 and options.laguerre

  for i in range(3):
    if info['size'][i] <= 0.0:                                                                        # any invalid size?
      info['size'][i] = float(info['grid'][i])/max(info['grid'])                                      # normalize to grid
      remarks.append('rescaling size {} to {}...'.format(['x','y','z'][i],info['size'][i]))

  if table.label_dimension(options.pos) != 3:
    errors.append('seed positions "{}" have dimension {}.'.format(options.pos,
                                                                  table.label_dimension(options.pos)))
  else:
    labels += [options.pos]

  if not options.normalized:               remarks.append('using real-space seed coordinates...')
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
  coords    = table.data[:,table.label_indexrange(options.pos)] * info['size'] if options.normalized \
        else  table.data[:,table.label_indexrange(options.pos)] - info['origin']
  eulers    = table.data[:,table.label_indexrange(options.eulers)] if hasEulers \
              else np.zeros(3*len(coords))
  grains    = table.data[:,table.label_indexrange(options.microstructure)].astype(int) if hasGrains \
              else 1+np.arange(len(coords))
  weights   = table.data[:,table.label_indexrange(options.weight)] if hasWeights \
              else np.zeros(len(coords))
  grainIDs  = np.unique(grains).astype('i')
  NgrainIDs = len(grainIDs)

# --- tessellate microstructure ------------------------------------------------------------

  x = (np.arange(info['grid'][0])+0.5)*info['size'][0]/info['grid'][0]
  y = (np.arange(info['grid'][1])+0.5)*info['size'][1]/info['grid'][1]
  z = (np.arange(info['grid'][2])+0.5)*info['size'][2]/info['grid'][2]
  
  damask.util.croak('tessellating...')

  grid = np.vstack(meshgrid2(x, y, z)).reshape(3,-1).T
  indices = laguerreTessellation(grid, coords, weights, grains, options.periodic, options.cpus)
    
  config_header = []
  formatwidth = 1+int(np.log10(NgrainIDs))

  if options.config:
    config_header += ['<microstructure>']
    for i,ID in enumerate(grainIDs):
      config_header += ['[Grain{}]'.format(str(ID).zfill(formatwidth)),
                        'crystallite 1',
                        '(constituent)\tphase {}\ttexture {}\tfraction 1.0'.format(options.phase,str(ID).rjust(formatwidth)),
                       ]
    if hasEulers:
      config_header += ['<texture>']
      theAxes = [] if options.axes is None else ['axes\t{} {} {}'.format(*options.axes)]
      for ID in grainIDs:
        eulerID = np.nonzero(grains == ID)[0][0]                                                    # find first occurrence of this grain id
        config_header += ['[Grain{}]'.format(str(ID).zfill(formatwidth)),
                          '(gauss)\tphi1 {:g}\tPhi {:g}\tphi2 {:g}'.format(*eulers[eulerID])
                         ] + theAxes
    config_header += ['<!skip>']
  
  header = [scriptID + ' ' + ' '.join(sys.argv[1:])] + config_header + ['origin x {} y {} z {}'.format(*info['origin'])]
  geom = damask.Geom(indices.reshape(info['grid'],order='F'),info['size'],options.homogenization,comments=header)
  damask.util.croak(geom)
  
  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(os.path.splitext(name)[0]+'.geom')
