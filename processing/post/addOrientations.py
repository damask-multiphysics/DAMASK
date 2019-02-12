#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os,sys
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add quaternion and/or Bunge Euler angle representation of crystal lattice orientation.
Orientation is given by quaternion, Euler angles, rotation matrix, or crystal frame coordinates
(i.e. component vectors of rotation matrix).
Additional (globally fixed) rotations of the lab frame and/or crystal frame can be applied.

""", version = scriptID)

representations = {
                  'quaternion': ['quat',4], #ToDo: Use here Rowenhorst names (qu/ro/om/ax?)
                  'rodrigues':  ['rodr',4],
                  'eulers':     ['eulr',3],
                  'matrix':     ['mtrx',9],
                  'angleaxis':  ['aaxs',4],
                  }

parser.add_option('-o',
                  '--output',
                  dest = 'output',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'output orientation formats {{{}}}'.format(', '.join(representations)))
parser.add_option('-d',
                  '--degrees',
                  dest = 'degrees',
                  action = 'store_true',
                  help = 'all angles in degrees')
parser.add_option('-R',
                  '--labrotation',
                  dest='labrotation',
                  type = 'float', nargs = 4, metavar = ' '.join(['float']*4),
                  help = 'angle and axis of additional lab frame rotation')
parser.add_option('-r',
                  '--crystalrotation',
                  dest='crystalrotation',
                  type = 'float', nargs = 4, metavar = ' '.join(['float']*4),
                  help = 'angle and axis of additional crystal frame rotation')
parser.add_option('--eulers',
                  dest = 'eulers',
                  metavar = 'string',
                  help = 'Euler angles label')
parser.add_option('--rodrigues',
                  dest = 'rodrigues',
                  metavar = 'string',
                  help = 'Rodrigues vector label')
parser.add_option('--matrix',
                  dest = 'matrix',
                  metavar = 'string',
                  help = 'orientation matrix label')
parser.add_option('--quaternion',
                  dest = 'quaternion',
                  metavar = 'string',
                  help = 'quaternion label')
parser.add_option('-x',
                  dest = 'x',
                  metavar = 'string',
                  help = 'label of lab x vector (expressed in crystal coords)')
parser.add_option('-y',
                  dest = 'y',
                  metavar = 'string',
                  help = 'label of lab y vector (expressed in crystal coords)')
parser.add_option('-z',
                  dest = 'z',
                  metavar = 'string',
                  help = 'label of lab z vector (expressed in crystal coords)')

parser.set_defaults(output = [],
                    labrotation     = (0.,1.,0.,0.),                                                # no rotation about 1,0,0
                    crystalrotation = (0.,1.,0.,0.),                                                # no rotation about 1,0,0
                    degrees = False,
                   )

(options, filenames) = parser.parse_args()

options.output = list(map(lambda x: x.lower(), options.output))
if options.output == [] or (not set(options.output).issubset(set(representations))):
  parser.error('output must be chosen from {}.'.format(', '.join(representations)))

input = [options.eulers     is not None,
         options.rodrigues  is not None,
         options.x          is not None and \
         options.y          is not None and \
         options.z          is not None,
         options.matrix     is not None,
         options.quaternion is not None,
        ]

if np.sum(input) != 1: parser.error('needs exactly one input format.')

(label,dim,inputtype) = [(options.eulers,representations['eulers'][1],'eulers'),
                         (options.rodrigues,representations['rodrigues'][1],'rodrigues'),
                         ([options.x,options.y,options.z],[3,3,3],'frame'),
                         (options.matrix,representations['matrix'][1],'matrix'),
                         (options.quaternion,representations['quaternion'][1],'quaternion'),
                        ][np.where(input)[0][0]]                                                    # select input label that was requested

crystalrotation = np.array(options.crystalrotation[1:4] + (options.crystalrotation[0],))            # Compatibility hack
labrotation     = np.array(options.labrotation[1:4],    + (options.labrotation[0],))                # Compatibility hack
r = damask.Rotation.fromAngleAxis(crystalrotation,options.degrees)                                  # crystal frame rotation
R = damask.Rotation.fromAngleAxis(labrotation,options.degrees)                                      #     lab frame rotation


# --- loop over input files ------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks -----------------------------------------

  errors  = []
  remarks = []
  
  if not np.all(table.label_dimension(label) == dim):  errors.append('input {} does not have dimension {}.'.format(label,dim))
  else:  column = table.label_index(label)

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  for output in options.output:
    if output in representations:
      table.labels_append(['{}_{}({})'.format(i+1,representations[output][0],label) \
                           for i in range(representations[output][1])])
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    if   inputtype == 'eulers':
      l = representations['eulers'][1]
      o = damask.Rotation.fromEulers(list(map(float,table.data[column:column+l])),options.degrees)
      
    elif inputtype == 'rodrigues':
      l = representations['rodrigues'][1]
      o = damask.Rotation.fromRodrigues(list(map(float,table.data[column:column+l])))
      
    elif inputtype == 'matrix':
      l = representations['matrix'][1]
      o = damask.Rotation.fromMatrix(list(map(float,table.data[column:column+l])))

    elif inputtype == 'frame':
      M = np.array(list(map(float,table.data[column[0]:column[0]+3] + \
                                  table.data[column[1]:column[1]+3] + \
                                  table.data[column[2]:column[2]+3]))).reshape(3,3).T
      o = damask.Rotation.fromMatrix(M/np.linalg.norm(M,axis=0))
      
    elif inputtype == 'quaternion':
      l = representations['quaternion'][1]
      o = damask.Rotation.fromQuaternion(list(map(float,table.data[column:column+l])))

    o= r*o*R                                                                                        # apply additional lab and crystal frame rotations

    for output in options.output:
      if   output == 'quaternion': table.data_append(o.asQuaternion())
      elif output == 'rodrigues':  table.data_append(o.asRodrigues())
      elif output == 'eulers':     table.data_append(o.asEulers(degrees=options.degrees))
      elif output == 'matrix':     table.data_append(o.asMatrix())
      elif output == 'angleaxis':  table.data_append(o.asAngleAxis(degrees=options.degrees))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
