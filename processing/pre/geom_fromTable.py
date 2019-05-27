#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser

import numpy as np

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Generate geometry description and material configuration from position, phase, and orientation (or microstructure) data.

""", version = scriptID)

parser.add_option('--coordinates',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'coordinates label (%default)')
parser.add_option('--phase',
                  dest = 'phase',
                  type = 'string', metavar = 'string',
                  help = 'phase label')
parser.add_option('--microstructure',
                  dest = 'microstructure',
                  type = 'string', metavar = 'string',
                  help = 'microstructure label')
parser.add_option('-q', '--quaternion',
                  dest = 'quaternion',
                  type = 'string', metavar='string',
                  help = 'quaternion label')
parser.add_option('--axes',
                  dest = 'axes',
                  type = 'string', nargs = 3, metavar = ' '.join(['string']*3),
                  help = 'orientation coordinate frame in terms of position coordinate frame [+x +y +z]')
parser.add_option('--homogenization',
                  dest = 'homogenization',
                  type = 'int', metavar = 'int',
                  help = 'homogenization index to be used [%default]')


parser.set_defaults(homogenization = 1,
                    pos            = 'pos',
                   )

(options,filenames) = parser.parse_args()

input = [ options.quaternion is not None,
         options.microstructure is not None,
        ]

if np.sum(input) != 1:
  parser.error('need either microstructure label or exactly one orientation input format.')
if options.axes is not None and not set(options.axes).issubset(set(['x','+x','-x','y','+y','-y','z','+z','-z'])):
  parser.error('invalid axes {} {} {}.'.format(*options.axes))

(label,dim,inputtype) = [(options.quaternion,4,'quaternion'),
                         (options.microstructure,1,'microstructure'),
                        ][np.where(input)[0][0]]                                                    # select input label that was requested


if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)
  table = damask.ASCIItable(name = name,readonly=True)
  table.head_read()                                                                                 # read ASCII header info

# ------------------------------------------ sanity checks ---------------------------------------  

  coordDim = table.label_dimension(options.pos)

  errors = []
  if not 3 >= coordDim >= 2:
    errors.append('coordinates "{}" need to have two or three dimensions.'.format(options.pos))
  if not np.all(table.label_dimension(label) == dim):
    errors.append('input "{}" needs to have dimension {}.'.format(label,dim))
  if options.phase and table.label_dimension(options.phase) != 1:
    errors.append('phase column "{}" is not scalar.'.format(options.phase))
  
  if errors != []:
    damask.util.croak(errors)
    continue

  table.data_readArray([options.pos] \
                       + (label if isinstance(label, list) else [label]) \
                       + ([options.phase] if options.phase else []))
  
  if coordDim == 2:
    table.data = np.insert(table.data,2,np.zeros(len(table.data)),axis=1)                           # add zero z coordinate for two-dimensional input
  if options.phase is None:
    table.data = np.column_stack((table.data,np.ones(len(table.data))))                             # add single phase if no phase column given

  grid,size = damask.util.coordGridAndSize(table.data[:,0:3])
  coords = [np.unique(table.data[:,i]) for i in range(3)]
  mincorner = np.array(list(map(min,coords)))
  origin = mincorner - 0.5*size/grid                                                                # shift from cell center to corner

  
# ------------------------------------------ process data ------------------------------------------

  colOri = table.label_index(label)+(3-coordDim)                                                    # column(s) of orientation data followed by 3 coordinates

  if inputtype == 'microstructure':

    grain = table.data[:,colOri]
    nGrains = len(np.unique(grain))

  elif inputtype == 'quaternion':

    colPhase = -1                                                                                   # column of phase data comes last
    index  = np.lexsort((table.data[:,0],table.data[:,1],table.data[:,2]))                          # index of position when sorting x fast, z slow
    grain = -np.ones(grid.prod(),dtype = int)                                                       # initialize empty microstructure
    orientations = []                                                                               # orientations
    multiplicity = []                                                                               # orientation multiplicity (number of group members)
    phases       = []                                                                               # phase info
    nGrains = 0                                                                                     # counter for detected grains
    existingGrains = np.arange(nGrains)
    myPos   = 0                                                                                     # position (in list) of current grid point


    for z in range(grid[2]):
      for y in range(grid[1]):
        for x in range(grid[0]):


          myData = table.data[index[myPos]]                                                         # read data for current grid point
          myPhase = int(myData[colPhase])
         
          o = damask.Rotation(myData[colOri:colOri+4])
          
          grain[myPos] = nGrains                                                                  # assign new grain to me ...
          nGrains += 1                                                                            # ... and update counter
          orientations.append(o)                                                                  # store new orientation for future comparison
          multiplicity.append(1)                                                                  # having single occurrence so far
          phases.append(myPhase)                                                                  # store phase info for future reporting
          existingGrains = np.arange(nGrains)                                                     # update list of existing grains

          myPos += 1

    
    grain += 1                                                                                    # offset from starting index 0 to 1
    

  formatwidth = 1+int(np.log10(nGrains))

  if inputtype == 'microstructure':
    config_header = []
  else:
    config_header = ['<microstructure>']
    for i,phase in enumerate(phases):
      config_header += ['[Grain%s]'%(str(i+1).zfill(formatwidth)),
                        'crystallite 1',
                        '(constituent)\tphase %i\ttexture %s\tfraction 1.0'%(phase,str(i+1).rjust(formatwidth)),
                       ]
  
    config_header += ['<texture>']
    for i,orientation in enumerate(orientations):
      config_header += ['[Grain%s]'%(str(i+1).zfill(formatwidth)),
                        '(gauss)\tphi1 %g\tPhi %g\tphi2 %g'%tuple(orientation.asEulers(degrees = True)),
                       ]
      if options.axes is not None: config_header += ['axes\t{} {} {}'.format(*options.axes)]
  
  header = [scriptID + ' ' + ' '.join(sys.argv[1:])] + config_header + ['origin x {} y {} z {}'.format(*origin)]
  geom = damask.Geom(grain.reshape(grid,order='F'),size,options.homogenization,comments=header)
  damask.util.croak(geom)
  
  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(os.path.splitext(name)[0]+'.geom')
