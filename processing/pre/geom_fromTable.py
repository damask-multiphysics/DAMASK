#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,re,string,math
import scipy.spatial, numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]



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
Generate geometry description and material configuration from position, phase, and orientation data.

""", version = scriptID)

parser.add_option('--coordinates',      dest='coordinates', type='string', metavar='string',
                  help='coordinates label')
parser.add_option('--phase',            dest='phase', type='string', metavar='string',
                  help='phase label')
parser.add_option('-t', '--tolerance', dest='tolerance', type='float', metavar='float',
                  help = 'angular tolerance for orientation squashing [%default]')
parser.add_option('-e', '--eulers',     dest='eulers', metavar='string',
                  help = 'Euler angles label')
parser.add_option('-d', '--degrees',    dest='degrees', action='store_true',
                  help = 'angles are given in degrees [%default]')
parser.add_option('-m', '--matrix',     dest='matrix', metavar='string',
                  help = 'orientation matrix label')
parser.add_option('-a',                 dest='a', metavar='string',
                  help = 'crystal frame a vector label')
parser.add_option('-b',                 dest='b', metavar='string',
                  help = 'crystal frame b vector label')
parser.add_option('-c',                 dest='c', metavar='string',
                  help = 'crystal frame c vector label')
parser.add_option('-q', '--quaternion', dest='quaternion', metavar='string',
                  help = 'quaternion label')
parser.add_option(      '--axes', dest='axes', nargs=3, metavar=' '.join(['string']*3),
                  help = 'orientation coordinate frame in terms of position coordinate frame [same]')
parser.add_option('-s', '--symmetry',   dest='symmetry', action='extend',
                  metavar='<string LIST>',
                  help = 'crystal symmetry [%s] {%s} '%(damask.Symmetry.lattices[-1],
                                                        ', '.join(damask.Symmetry.lattices[1:])))
parser.add_option('--homogenization',   dest='homogenization', type='int', metavar='int',
                  help='homogenization index to be used [%default]')
parser.add_option('--crystallite',      dest='crystallite', type='int', metavar='int',
                  help='crystallite index to be used [%default]')

parser.set_defaults(symmetry       = [damask.Symmetry.lattices[-1]],
                    tolerance      = 0.0,
                    degrees        = False,
                    homogenization = 1,
                    crystallite    = 1,
                   )

(options,filenames) = parser.parse_args()

input = [options.eulers     != None,
         options.a          != None and \
         options.b          != None and \
         options.c          != None,
         options.matrix     != None,
         options.quaternion != None,
        ]

if np.sum(input) != 1: parser.error('needs exactly one orientation input format...')
if options.axes != None and not set(options.axes).issubset(set(['x','+x','-x','y','+y','-y','z','+z','-z'])):
  parser.error('invalid axes %s %s %s'%tuple(options.axes))

(label,dim,inputtype) = [(options.eulers,3,'eulers'),
                         ([options.a,options.b,options.c],[3,3,3],'frame'),
                         (options.matrix,9,'matrix'),
                         (options.quaternion,4,'quaternion'),
                        ][np.where(input)[0][0]]                                                    # select input label that was requested
toRadians = math.pi/180.0 if options.degrees else 1.0                                               # rescale degrees to radians
options.tolerance *= toRadians                                                                      # angular tolerance in radians

# --- loop over input files -------------------------------------------------------------------------
if filenames == []:
  filenames = ['STDIN']

for name in filenames:
  if name == 'STDIN':
    file = {'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m\n')
  else:
    if not os.path.exists(name): continue
    file = {'name':name,
            'input':open(name),
            'output':open(name + '_tmp','w'),
            'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],buffered=False)                            # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info

  errors = []
  if table.label_dimension(options.coordinates) != 2:
    errors.append('coordinates %s need to have two dimensions...'%options.coordinates)
  if not np.all(table.label_dimension(label) == dim):
    errors.append('orientation %s needs to have dimension %i...\n'%(label,dim))
  if options.phase != None and table.label_dimension(options.phase) != 1:
    errors.append('phase column %s is not scalar...'%options.phase)
  
  if errors == []:
    table.data_readArray([options.coordinates,label]+([] if options.phase == None else [options.phase]))
    
    if options.phase == None:
      table.data = np.column_stack((table.data,np.ones(len(table.data))))                         # add single phase if no phase column given
    
    coordsX = np.unique(table.data[:,0])
    coordsY = np.unique(table.data[:,1])
    nX = len(coordsX)
    nY = len(coordsY)
    dX = (coordsX[-1]-coordsX[0])/(nX-1)
    dY = (coordsY[-1]-coordsY[0])/(nY-1)

    if    nX*nY != len(table.data) \
       or np.any(np.abs(np.log10((coordsX[1:]-coordsX[:-1])/dX)) > 0.01) \
       or np.any(np.abs(np.log10((coordsY[1:]-coordsY[:-1])/dY)) > 0.01):
      errors.append('data is not on square grid...')
  
  if errors != []:
    file['croak'].write('\n'.join(errors)+'\n')
    table.close(dismiss = True)
    continue
  
# ------------------------------------------ process data ------------------------------------------
  
  colOri = table.label_index(label)                                                               # column(s) of orientation data
  colPhase = colOri + np.sum(dim)                                                                 # column of phase data comes after orientation
  index = np.lexsort((table.data[:,0],table.data[:,1]))                                           # index of rank when sorting x fast, y slow
  rank  = np.argsort(index)                                                                       # rank of index
  KDTree = scipy.spatial.KDTree((table.data[:,:2]-np.array([coordsX[0],coordsY[0]])) \
                                /               np.array([dX,dY]))                                # build KDTree with dX = dY = 1
  
  microstructure = np.zeros(nX*nY,dtype='uint32')                                                 # initialize empty microstructure
  symQuats = []                                                                                   # empty list of sym equiv orientations
  phases   = []                                                                                   # empty list of phase info
  nGrains = 0                                                                                     # counter for detected grains
  myRank  = 0                                                                                     # rank of current grid point
  for y in xrange(nY):
    for x in xrange(nX):
      if (myRank+1)%max(1,nX*nY/100) == 0:
        file['croak'].write('.')
      myData = table.data[index[myRank]]
      mySym = options.symmetry[min(int(myData[colPhase]),len(options.symmetry))-1]                # select symmetry from option (take last specified option for all with higher index)
      if inputtype == 'eulers':
        o = damask.Orientation(Eulers=toRadians*\
                               np.array(map(float,myData[colOri:colOri+3])),
                               symmetry=mySym).reduced()
      elif inputtype == 'matrix':
        o = damask.Orientation(matrix=\
                               np.array([map(float,myData[colOri:colOri+9])]).reshape(3,3).transpose(),
                               symmetry=mySym).reduced()
      elif inputtype == 'frame':
        o = damask.Orientation(matrix=\
                               np.array([map(float,myData[colOri[0]:colOri[0]+3] + \
                                                   myData[colOri[1]:colOri[1]+3] + \
                                                   myData[colOri[2]:colOri[2]+3]
                                       )]).reshape(3,3),
                               symmetry=mySym).reduced()
      elif inputtype == 'quaternion':
        o = damask.Orientation(quaternion=\
                               np.array(map(float,myData[colOri:colOri+4])),
                               symmetry=mySym).reduced()

      oInv = o.quaternion.conjugated()
      neighbors = KDTree.query_ball_point([x,y], 3)                                               # search points within radius
      breaker = False

      for n in neighbors:                                                                         # check each neighbor
        if myRank <= rank[n] or table.data[n,colPhase] != myData[colPhase]: continue              # skip myself, anyone further ahead (cannot yet have a grain ID), and other phases
        for q in symQuats[microstructure[rank[n]]-1]:
          if abs((q*oInv).asAngleAxis()[0]) <= options.tolerance:                                      # found existing orientation resembling me
            microstructure[myRank] = microstructure[rank[n]]
            breaker = True; break
        if breaker: break

      if microstructure[myRank] == 0:                                                             # no other orientation resembled me
        nGrains += 1
        microstructure[myRank] = nGrains
        symQuats.append(o.equivalentQuaternions())                                                # store all symmetrically equivalent orientations for future comparison
        phases.append(myData[colPhase])                                                           # store phase info for future reporting

      myRank += 1

  file['croak'].write('\n')
#--- generate header ----------------------------------------------------------------------------

  info = {
          'grid':    np.array([nX,nY,1]),
          'size':    np.array([coordsX[-1]-coordsX[0],
                               coordsY[-1]-coordsY[0],
                               min((coordsX[-1]-coordsX[0])/nX,
                                   (coordsY[-1]-coordsY[0])/nY,
                                  )
                              ]),
          'origin':  np.array([coordsX[0],coordsY[0],0.0]),
          'microstructures': nGrains,
          'homogenization':  options.homogenization,
         }
  
  file['croak'].write('grid     a b c: %s\n'%(' x '.join(map(str,info['grid']))) + \
                      'size     x y z: %s\n'%(' x '.join(map(str,info['size']))) + \
                      'origin   x y z: %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization: %i\n'%info['homogenization'] + \
                      'microstructures: %i\n'%info['microstructures'])
  
#--- write header ---------------------------------------------------------------------------------

  formatwidth = 1+int(math.log10(info['microstructures']))

  config_header = ['<microstructure>']
  for i,phase in enumerate(phases):
    config_header += ['[Grain%s]'%(str(i+1).zfill(formatwidth)),
                      'crystallite %i'%options.crystallite,
                      '(constituent)\tphase %i\ttexture %s\tfraction 1.0'%(phase,str(i+1).rjust(formatwidth)),
                     ]
  
  config_header += ['<texture>']
  for i,quats in enumerate(symQuats):
    config_header += ['[Grain%s]'%(str(i+1).zfill(formatwidth)),
                      'axes\t%s %s %s'%tuple(options.axes) if options.axes != None else '',
                      '(gauss)\tphi1 %g\tPhi %g\tphi2 %g\tscatter 0.0\tfraction 1.0'%tuple(np.degrees(quats[0].asEulers())),
                     ]
  
  table.labels_clear()
  table.info_clear()
  table.info_append([
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
    "size\tx %f\ty %f\tz %f"%(info['size'][0],info['size'][1],info['size'][2],),
    "origin\tx %f\ty %f\tz %f"%(info['origin'][0],info['origin'][1],info['origin'][2],),
    "homogenization\t%i"%info['homogenization'],
    "microstructures\t%i"%(info['microstructures']),
    config_header,
    ])
  table.head_write()
  
# --- write microstructure information ------------------------------------------------------------

  table.data = microstructure.reshape(info['grid'][1]*info['grid'][2],info['grid'][0])
  table.data_writeArray('%%%ii'%(formatwidth),delimiter=' ')
  
#--- output finalization --------------------------------------------------------------------------

  table.close()
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',
              os.path.splitext(file['name'])[0] + '.geom')
