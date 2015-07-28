#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,string
import numpy as np
from optparse import OptionParser
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

def laguerreTessellation(undeformed, coords, weights):

    weight = np.power(np.tile(weights, 27),2)                        # Laguerre weights (squared)
    micro = np.zeros(undeformed.shape[0])
    N = coords.shape[0]                                              # Number of seeds points 
    periodic =  np.array([
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
                         ]).astype(float)
    
    for i,vec in enumerate(periodic):
        seeds = np.append(seeds, coords+vec, axis=0) if i > 0 else coords+vec
    
    for i,point in enumerate(undeformed):
        
        tmp = np.repeat(point.reshape(3,1), N*27, axis=1).T
        dist = np.sum((tmp - seeds)*(tmp - seeds),axis=1) - weight
        micro[i] = np.argmin(dist)%N + 1
    
    return micro

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
identifiers = {
        'grid':   ['a','b','c'],
        'size':   ['x','y','z'],
        'origin': ['x','y','z'],
          }
mappings = {
        'grid':            lambda x: int(x),
        'size':            lambda x: float(x),
        'origin':          lambda x: float(x),
        'homogenization':  lambda x: int(x),
        'microstructures': lambda x: int(x),
          }

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Generate geometry description and material configuration by standard Voronoi tessellation of given seeds file.

""", version = scriptID)

parser.add_option('-g', '--grid', dest='grid', type='int', nargs = 3, metavar=' '.join(['int']*3),
                  help='a,b,c grid of hexahedral box [from seeds file]')
parser.add_option('-s', '--size', dest='size', type='float', nargs = 3, metavar=' '.join(['float']*3),
                  help='x,y,z size of hexahedral box [1.0 along largest grid point number]')
parser.add_option('-o', '--origin', dest='origin', type='float', nargs = 3, metavar=' '.join(['float']*3),
                  help='offset from old to new origin of grid')
parser.add_option('--homogenization', dest='homogenization', type='int', metavar = 'int',
                  help='homogenization index to be used [%default]')
parser.add_option('--phase', dest='phase', type='int', metavar = 'int',
                  help='phase index to be used [%default]')
parser.add_option('--crystallite', dest='crystallite', type='int', metavar = 'int',
                  help='crystallite index to be used [%default]')
parser.add_option('-c', '--configuration', dest='config', action='store_true',
                  help='output material configuration [%default]')
parser.add_option('-r', '--rnd', dest='randomSeed', type='int', metavar='int',
                  help='seed of random number generator for second phase distribution [%default]')
parser.add_option('--secondphase', type='float', dest='secondphase', metavar= 'float',
                  help='volume fraction of randomly distribute second phase [%default]')
parser.add_option('-l', '--laguerre', dest='laguerre', action='store_true',
                  help='use Laguerre (weighted Voronoi) tessellation [%default]')
parser.set_defaults(grid   = (0,0,0),
                    size   = (0.0,0.0,0.0),
                    origin = (0.0,0.0,0.0),
                    homogenization = 1,
                    phase          = 1,
                    crystallite    = 1,
                    secondphase    = 0.0,
                    config = False,
                    laguerre = False,
                    randomSeed = None,
                  )
(options,filenames) = parser.parse_args()

if options.secondphase > 1.0 or options.secondphase < 0.0:
 parser.error('volume fraction of second phase (%f) out of bounds...'%options.secondphase)

# --- loop over input files -------------------------------------------------------------------------
if filenames == []:
  filenames = ['STDIN']

for name in filenames:
  if name == 'STDIN':
    file = {'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m\n')
  else:
    if not os.path.exists(name): continue
    file = {'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],buffered=False)                            # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info

  labels = []
  if np.all(table.label_index(['1_coords','2_coords','3_coords']) != -1):
    coords = ['1_coords','2_coords','3_coords']
  elif np.all(table.label_index(['x','y','z']) != -1):
    coords = ['x','y','z']
  else:
    file['croak'].write('no coordinate data (1/2/3_coords | x/y/z) found ...')
    continue

  labels += coords
  hasEulers = np.all(table.label_index(['phi1','Phi','phi2']) != -1)
  if hasEulers:
    labels += ['phi1','Phi','phi2']
    
  hasGrains = table.label_index('microstructure') != -1
  if hasGrains:
    labels += ['microstructure']
    
  hasWeight = table.label_index('weight') != -1
  if hasWeight:
    labels += ['weight']
      
  table.data_readArray(labels)
  coords = table.data[:,table.label_index(coords)]
  eulers = table.data[:,table.label_index(['phi1','Phi','phi2'])] if hasEulers else np.zeros(3*len(coords))
  grain = table.data[:,table.label_index('microstructure')]       if hasGrains else 1+np.arange(len(coords))
  weights = table.data[:,table.label_index('weight')]             if hasWeight else np.zeros(len(coords))
  grainIDs = np.unique(grain).astype('i')


#--- interpret header ----------------------------------------------------------------------------
  info = {
          'grid':    np.zeros(3,'i'),
          'size':    np.array(options.size),
          'origin':  np.zeros(3,'d'),
          'microstructures':  0,
          'homogenization': options.homogenization,
         }
  newInfo = {
          'microstructures': 0,
         }
  extra_header = []

  for header in table.info:
    headitems = map(str.lower,header.split())
    if len(headitems) == 0: continue
    if headitems[0] in mappings.keys():
      if headitems[0] in identifiers.keys():
        for i in xrange(len(identifiers[headitems[0]])):
          info[headitems[0]][i] = \
            mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])
      else:
        info[headitems[0]] = mappings[headitems[0]](headitems[1])
    else:
      extra_header.append(header)

  if info['microstructures'] != len(grainIDs):
    file['croak'].write('grain data not matching grain count (%i)...\n'%(len(grainIDs)))
    info['microstructures'] = len(grainIDs)
  
  if 0 not in options.grid:                                                                         # user-specified grid
    info['grid'] = np.array(options.grid)

  for i in xrange(3):
    if info['size'][i] <= 0.0:                                                                      # any invalid size?
      info['size'][i] = float(info['grid'][i])/max(info['grid'])
      file['croak'].write('rescaling size %s...\n'%{0:'x',1:'y',2:'z'}[i])

  file['croak'].write('grains to map:  %i\n'%info['microstructures'] + \
                      'grid     a b c: %s\n'%(' x '.join(map(str,info['grid']))) + \
                      'size     x y z: %s\n'%(' x '.join(map(str,info['size']))) + \
                      'origin   x y z: %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization: %i\n'%info['homogenization'])
  
  if np.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    continue
  if np.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    continue
  if info['microstructures'] == 0:
    file['croak'].write('no grain info found.\n')
    continue

#--- prepare data ---------------------------------------------------------------------------------
  eulers = eulers.T

#--- switch according to task ---------------------------------------------------------------------
  if options.config:                                                                                # write config file
    phase = np.empty(info['microstructures'],'i')
    phase.fill(options.phase)
    phase[0:int(float(info['microstructures'])*options.secondphase)] = options.phase+1
    randomSeed = int(os.urandom(4).encode('hex'), 16)  if options.randomSeed == None else options.randomSeed         # radom seed per file for second phase
    np.random.seed(randomSeed)
    np.random.shuffle(phase)
    formatwidth = 1+int(math.log10(info['microstructures']))
    file['output'].write('#' + scriptID + ' ' + ' '.join(sys.argv[1:])+'\n')
    if options.secondphase > 0.0: file['output'].write('# random seed for second phase %i\n'%randomSeed)
    file['output'].write('\n<microstructure>\n')
    for i,ID in enumerate(grainIDs):
      file['output'].write('\n[Grain%s]\n'%(str(ID).zfill(formatwidth)) + \
                           'crystallite %i\n'%options.crystallite + \
                           '(constituent)\tphase %i\ttexture %s\tfraction 1.0\n'%(phase[i],str(ID).rjust(formatwidth)))
  
    file['output'].write('\n<texture>\n')
    for ID in grainIDs:
      eulerID = np.nonzero(grain == ID)[0][0]                                                       # find first occurrence of this grain id
      file['output'].write('\n[Grain%s]\n'%(str(ID).zfill(formatwidth)) + \
                           '(gauss)\tphi1 %g\tPhi %g\tphi2 %g\tscatter 0.0\tfraction 1.0\n'%(eulers[0,eulerID],
                                                                                             eulers[1,eulerID],
                                                                                             eulers[2,eulerID]))

  else:                                                                                             # write geometry file
    x = (np.arange(info['grid'][0])+0.5)*info['size'][0]/info['grid'][0]
    y = (np.arange(info['grid'][1])+0.5)*info['size'][1]/info['grid'][1]
    z = (np.arange(info['grid'][2])+0.5)*info['size'][2]/info['grid'][2]
    
    if not options.laguerre:
      coords = (coords*info['size']).T
      undeformed = np.vstack(map(np.ravel, meshgrid2(x, y, z)))
  
      file['croak'].write('tessellating...\n')
      indices = damask.core.math.periodicNearestNeighbor(\
                info['size'],\
                np.eye(3),\
                undeformed,coords)//3**3 + 1                                                          # floor division to kill periodic images
      indices = grain[indices-1]
    else :
      undeformed = np.vstack(np.meshgrid(x, y, z)).reshape(3,-1).T
      indices = laguerreTessellation(undeformed, coords, weights)

    newInfo['microstructures'] = info['microstructures']
    for i in grainIDs:
      if i not in indices: newInfo['microstructures'] -= 1
    file['croak'].write(('all' if newInfo['microstructures']  == info['microstructures'] else 'only') +
                        ' %i'%newInfo['microstructures'] + 
                        ('' if newInfo['microstructures'] == info['microstructures'] else \
                        ' out of %i'%info['microstructures']) +
                        ' grains mapped.\n')

#--- write header ---------------------------------------------------------------------------------
    table.labels_clear()
    table.info_clear()
    table.info_append(extra_header+[
      scriptID + ' ' + ' '.join(sys.argv[1:]),
      "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
      "size\tx %f\ty %f\tz %f"%(info['size'][0],info['size'][1],info['size'][2],),
      "origin\tx %f\ty %f\tz %f"%(info['origin'][0],info['origin'][1],info['origin'][2],),
      "homogenization\t%i"%info['homogenization'],
      "microstructures\t%i"%(newInfo['microstructures']),
      ])
    table.head_write()
    
# --- write microstructure information ------------------------------------------------------------
    formatwidth = 1+int(math.log10(newInfo['microstructures']))
    table.data = indices.reshape(info['grid'][1]*info['grid'][2],info['grid'][0])
    table.data_writeArray('%%%ii'%(formatwidth),delimiter=' ')
    
#--- output finalization --------------------------------------------------------------------------

  table.close()  
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',
              os.path.splitext(file['name'])[0] +'%s'%('_material.config' if options.config else '.geom'))
