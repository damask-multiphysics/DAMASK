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
Converts ASCII table. Input can be microstructure or orientation (as quaternion). For the latter,
phase information can be given additionally.

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

input = [options.quaternion     is not None,
         options.microstructure is not None,
        ]

if np.sum(input) != 1:
  parser.error('need either microstructure or quaternion (and optionally phase) as input.')
if options.microstructure is not None and options.phase is not None:
  parser.error('need either microstructure or phase (and mandatory quaternion) as input.')
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


  indices  = np.lexsort((table.data[:,0],table.data[:,1],table.data[:,2]))                          # indices of position when sorting x fast, z slow
  microstructure = np.empty(grid,dtype = int)                                                       # initialize empty microstructure
  i = 0

  if inputtype == 'microstructure':
    for z in range(grid[2]):
      for y in range(grid[1]):
        for x in range(grid[0]):
          microstructure[x,y,z] = table.data[indices[i],3]
          i+=1

    config_header = []

  elif inputtype == 'quaternion':
    unique,unique_inverse = np.unique(table.data[:,3:8],return_inverse=True,axis=0)

    for z in range(grid[2]):
      for y in range(grid[1]):
        for x in range(grid[0]):
          microstructure[x,y,z] = unique_inverse[indices[i]]+1
          i+=1

    config_header = ['<texture>']
    for i,data in enumerate(unique):
      ori = damask.Rotation(data[0:4])
      config_header += ['[Grain{}]'.format(i+1),
                        '(gauss)\tphi1 {:g}\tPhi {:g}\tphi2 {:g}'.format(*ori.asEulers(degrees = True)),
                       ]
      if options.axes is not None: config_header += ['axes\t{} {} {}'.format(*options.axes)]

    config_header += ['<microstructure>']
    for i,data in enumerate(unique):
      config_header += ['[Grain{}]'.format(i+1),
                        'crystallite 1',
                        '(constituent)\tphase {}\ttexture {}\tfraction 1.0'.format(int(data[4]),i+1),
                       ]

  header = [scriptID + ' ' + ' '.join(sys.argv[1:])]\
         + config_header
  geom = damask.Geom(microstructure,size,origin,
                     homogenization=options.homogenization,comments=header)
  damask.util.croak(geom)

  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(os.path.splitext(name)[0]+'.geom')
