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

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Generate geometry description and material configuration from position, phase, and orientation (or microstructure) data.

""", version = scriptID)

parser.add_option('--coordinates',
                  dest = 'coordinates',
                  type = 'string', metavar = 'string',
                  help = 'coordinates label')
parser.add_option('--phase',
                  dest = 'phase',
                  type = 'string', metavar = 'string',
                  help = 'phase label')
parser.add_option('--microstructure',
                  dest = 'microstructure',
                  type = 'string', metavar = 'string',
                  help = 'microstructure label')
parser.add_option('-t', '--tolerance',
                  dest = 'tolerance',
                  type = 'float', metavar = 'float',
                  help = 'angular tolerance for orientation squashing [%default]')
parser.add_option('-e', '--eulers',
                  dest = 'eulers',
                  type = 'string', metavar = 'string',
                  help = 'Euler angles label')
parser.add_option('-d', '--degrees',
                  dest = 'degrees',
                  action = 'store_true',
                  help = 'angles are given in degrees [%default]')
parser.add_option('-m', '--matrix',
                  dest = 'matrix',
                  type = 'string', metavar = 'string',
                  help = 'orientation matrix label')
parser.add_option('-a',
                  dest='a',
                  type = 'string', metavar = 'string',
                  help = 'crystal frame a vector label')
parser.add_option('-b',
                  dest='b',
                  type = 'string', metavar = 'string',
                  help = 'crystal frame b vector label')
parser.add_option('-c',
                  dest = 'c',
                  type =  'string', metavar='string',
                  help = 'crystal frame c vector label')
parser.add_option('-q', '--quaternion',
                  dest = 'quaternion',
                  type = 'string', metavar='string',
                  help = 'quaternion label')
parser.add_option('--axes',
                  dest = 'axes',
                  type = 'string', nargs = 3, metavar = ' '.join(['string']*3),
                  help = 'orientation coordinate frame in terms of position coordinate frame [same]')
parser.add_option('-s', '--symmetry',
                  dest = 'symmetry',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'crystal symmetry %default {{{}}} '.format(', '.join(damask.Symmetry.lattices[1:])))
parser.add_option('--homogenization',
                  dest = 'homogenization',
                  type = 'int', metavar = 'int',
                  help = 'homogenization index to be used [%default]')
parser.add_option('--crystallite',
                  dest = 'crystallite',
                  type = 'int', metavar = 'int',
                  help = 'crystallite index to be used [%default]')

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

if np.sum(input) != 1 and options.microstructure == None:
  parser.error('need either microstructure label or exactly one orientation input format.')
if options.axes != None and not set(options.axes).issubset(set(['x','+x','-x','y','+y','-y','z','+z','-z'])):
  parser.error('invalid axes {} {} {}.'.format(*options.axes))

(label,dim,inputtype) = [(options.eulers,3,'eulers'),
                         ([options.a,options.b,options.c],[3,3,3],'frame'),
                         (options.matrix,9,'matrix'),
                         (options.quaternion,4,'quaternion'),
                        ][np.where(input)[0][0]]                                                    # select input label that was requested
toRadians = math.pi/180.0 if options.degrees else 1.0                                               # rescale degrees to radians
options.tolerance *= toRadians                                                                      # ensure angular tolerance in radians

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name, outname = os.path.splitext(name)[0] + '.geom',
                              buffered = False)
  except: continue
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name else ''))

# ------------------------------------------ read head ---------------------------------------  

  table.head_read()                                                                                 # read ASCII header info

# ------------------------------------------ sanity checks ---------------------------------------  

  errors = []
  if not 3 >= table.label_dimension(options.coordinates) >= 2:                                      # TODO need to deal with 3D case!!
    errors.append('coordinates {} need to have two or three dimensions.'.format(options.coordinates))
  if not np.all(table.label_dimension(label) == dim):
    errors.append('orientation {} needs to have dimension {}.'.format(label,dim))
  if options.phase != None and table.label_dimension(options.phase) != 1:
    errors.append('phase column {} is not scalar.'.format(options.phase))
  
  if errors == []:                                                                                # so far no errors?
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
      errors.append('data is not on square grid.')
  
  if errors != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue
  
# ------------------------------------------ process data ------------------------------------------
  
  colOri = table.label_index(label)                                                               # column(s) of orientation data
  colPhase = colOri + np.sum(dim)                                                                 # column of phase data comes after orientation
  index = np.lexsort((table.data[:,0],table.data[:,1]))                                           # index of rank when sorting x fast, y slow
  rank  = np.argsort(index)                                                                       # rank of index
  KDTree = scipy.spatial.KDTree((table.data[:,:2]-np.array([coordsX[0],coordsY[0]])) \
                                /              np.array([dX,dY]))                                 # build KDTree with dX = dY = 1
  
  microstructure = np.zeros(nX*nY,dtype = 'uint32')                                               # initialize empty microstructure
  symQuats = []                                                                                   # empty list of sym equiv orientations
  phases   = []                                                                                   # empty list of phase info
  nGrains = 0                                                                                     # counter for detected grains
  myRank  = 0                                                                                     # rank of current grid point
  for y in xrange(nY):
    for x in xrange(nX):
      if (myRank+1)%(nX*nY/100.) < 1: table.croak('.',False)
      myData = table.data[index[myRank]]
      mySym = options.symmetry[min(int(myData[colPhase]),len(options.symmetry))-1]                # select symmetry from option (take last specified option for all with higher index)
      if inputtype == 'eulers':
        o = damask.Orientation(Eulers = np.array(map(float,myData[colOri:colOri+3]))*toRadians,
                               symmetry = mySym).reduced()
      elif inputtype == 'matrix':
        o = damask.Orientation(matrix = np.array([map(float,myData[colOri:colOri+9])]).reshape(3,3).transpose(),
                               symmetry = mySym).reduced()
      elif inputtype == 'frame':
        o = damask.Orientation(matrix = np.array([map(float,myData[colOri[0]:colOri[0]+3] + \
                                                            myData[colOri[1]:colOri[1]+3] + \
                                                            myData[colOri[2]:colOri[2]+3]
                                                     )]).reshape(3,3),
                               symmetry = mySym).reduced()
      elif inputtype == 'quaternion':
        o = damask.Orientation(quaternion = np.array(map(float,myData[colOri:colOri+4])),
                               symmetry = mySym).reduced()

      oInv = o.quaternion.conjugated()
      neighbors = KDTree.query_ball_point([x,y], 3)                                               # search points within radius
      breaker = False

      for n in neighbors:                                                                         # check each neighbor
        if myRank <= rank[n] or table.data[n,colPhase] != myData[colPhase]: continue              # skip myself, anyone further ahead (cannot yet have a grain ID), and other phases
        for q in symQuats[microstructure[rank[n]]-1]:
          if abs((q*oInv).asAngleAxis()[0]) <= options.tolerance:                                 # found existing orientation resembling me
            microstructure[myRank] = microstructure[rank[n]]
            breaker = True; break
        if breaker: break

      if microstructure[myRank] == 0:                                                             # no other orientation resembled me
        nGrains += 1                                                                              # make new grain ...
        microstructure[myRank] = nGrains                                                          # ... and assign to me
        symQuats.append(o.equivalentQuaternions())                                                # store all symmetrically equivalent orientations for future comparison
        phases.append(myData[colPhase])                                                           # store phase info for future reporting

      myRank += 1

  table.croak('')

# --- generate header ----------------------------------------------------------------------------

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

  table.croak(['grid     a b c:  %s'%(' x '.join(map(str,info['grid']))),
               'size     x y z:  %s'%(' x '.join(map(str,info['size']))),
               'origin   x y z:  %s'%(' : '.join(map(str,info['origin']))),
               'homogenization:  %i'%info['homogenization'],
               'microstructures: %i'%info['microstructures'],
              ])
    
# --- write header ---------------------------------------------------------------------------------

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
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=info['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=info['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=info['origin']),
    "homogenization\t{homog}".format(homog=info['homogenization']),
    "microstructures\t{microstructures}".format(microstructures=info['microstructures']),
    config_header,
    ])
  table.head_write()
  
# --- write microstructure information ------------------------------------------------------------

  table.data = microstructure.reshape(info['grid'][1]*info['grid'][2],info['grid'][0])
  table.data_writeArray('%%%ii'%(formatwidth),delimiter=' ')
  
#--- output finalization --------------------------------------------------------------------------

  table.close()
