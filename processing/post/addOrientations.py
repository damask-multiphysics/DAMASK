#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                     convention conformity checks
# --------------------------------------------------------------------

def check_Eulers(eulers):
  if np.any(eulers < 0.0) or np.any(eulers > 2.0*math.pi) or eulers[1] > math.pi:                   # Euler angles within valid range?
    raise ValueError('Euler angles outside of [0..2π],[0..π],[0..2π].\n{} {} {}.'.format(*eulers))
  return eulers

def check_quaternion(q):
  if q[0] < 0.0:                                                                                    # positive first quaternion component?
    raise ValueError('quaternion has negative first component.\n{}'.format(q))
  if not(np.isclose(np.linalg.norm(q), 1.0)):                                                       # unit quaternion?
    raise ValueError('quaternion is not of unit length.\n{} {} {} {}'.format(*q))
  return q
 
def check_matrix(M):
  if abs(1.0-np.linalg.det(M)) > 1e-8:                                                              # proper rotation?
    raise ValueError('matrix is not a proper rotation.\n{}'.format(M))
  if    abs(np.dot(M[0],M[1]))    > 1e-8 \
     or abs(np.dot(M[1],M[2]))    > 1e-8 \
     or abs(np.dot(M[2],M[0]))    > 1e-8:                                                           # all orthogonal?
    raise ValueError('matrix is not orthogonal.\n{}'.format(M))
  return M

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add quaternion and/or Bunge Euler angle representation of crystal lattice orientation.
Orientation is given by quaternion, Euler angles, rotation matrix, or crystal frame coordinates
(i.e. component vectors of rotation matrix).
Additional (globally fixed) rotations of the lab frame and/or crystal frame can be applied.

""", version = scriptID)

outputChoices = {
                  'quaternion': ['quat',4],
                  'rodrigues':  ['rodr',3],
                  'eulers':     ['eulr',3],
                  'matrix':     ['mtrx',9],
                  'angleaxis':  ['aaxs',4],
                }

parser.add_option('-o', '--output',
                  dest = 'output',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'output orientation formats {{{}}}'.format(', '.join(outputChoices)))
parser.add_option('-d', '--degrees',
                  dest = 'degrees',
                  action = 'store_true',
                  help = 'all angles in degrees')
parser.add_option('-R', '--labrotation',
                  dest='labrotation',
                  type = 'float', nargs = 4, metavar = ' '.join(['float']*4),
                  help = 'angle and axis of additional lab frame rotation')
parser.add_option('-r', '--crystalrotation',
                  dest='crystalrotation',
                  type = 'float', nargs = 4, metavar = ' '.join(['float']*4),
                  help = 'angle and axis of additional crystal frame rotation')
parser.add_option(      '--eulers',
                  dest = 'eulers',
                  type = 'string', metavar = 'string',
                  help = 'Euler angles label')
parser.add_option(      '--rodrigues',
                  dest = 'rodrigues',
                  type = 'string', metavar = 'string',
                  help = 'Rodrigues vector label')
parser.add_option(      '--matrix',
                  dest = 'matrix',
                  type = 'string', metavar = 'string',
                  help = 'orientation matrix label')
parser.add_option(       '--quaternion',
                  dest = 'quaternion',
                  type = 'string', metavar = 'string',
                  help = 'quaternion label')
parser.add_option('-x',
                  dest = 'x',
                  type = 'string', metavar = 'string',
                  help = 'label of lab x vector (expressed in crystal coords)')
parser.add_option('-y',
                  dest = 'y',
                  type = 'string', metavar = 'string',
                  help = 'label of lab y vector (expressed in crystal coords)')
parser.add_option('-z',
                  dest = 'z',
                  type = 'string', metavar = 'string',
                  help = 'label of lab z vector (expressed in crystal coords)')

parser.set_defaults(output = [],
                    labrotation     = (0.,1.,1.,1.),                                                # no rotation about 1,1,1
                    crystalrotation = (0.,1.,1.,1.),                                                # no rotation about 1,1,1
                    degrees = False,
                   )

(options, filenames) = parser.parse_args()

options.output = list(map(lambda x: x.lower(), options.output))
if options.output == [] or (not set(options.output).issubset(set(outputChoices))):
  parser.error('output must be chosen from {}.'.format(', '.join(outputChoices)))

input = [options.eulers     is not None,
         options.rodrigues  is not None,
         options.x          is not None and \
         options.y          is not None and \
         options.z          is not None,
         options.matrix     is not None,
         options.quaternion is not None,
        ]

if np.sum(input) != 1: parser.error('needs exactly one input format.')

(label,dim,inputtype) = [(options.eulers,3,'eulers'),
                         (options.rodrigues,3,'rodrigues'),
                         ([options.x,options.y,options.z],[3,3,3],'frame'),
                         (options.matrix,9,'matrix'),
                         (options.quaternion,4,'quaternion'),
                        ][np.where(input)[0][0]]                                                    # select input label that was requested

toRadians = math.pi/180.0 if options.degrees else 1.0                                               # rescale degrees to radians
r = damask.Quaternion.fromAngleAxis(toRadians*options.crystalrotation[0],options.crystalrotation[1:]) # crystal frame rotation
R = damask.Quaternion.fromAngleAxis(toRadians*options.    labrotation[0],options.    labrotation[1:]) #     lab frame rotation

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
    if output in outputChoices:
      table.labels_append(['{}_{}({})'.format(i+1,outputChoices[output][0],label) \
                           for i in range(outputChoices[output][1])])
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    if   inputtype == 'eulers':
      
      o = damask.Orientation(Eulers = check_Eulers(np.array(list(map(float,table.data[column:column+3])))*toRadians))
    elif inputtype == 'rodrigues':
      o = damask.Orientation(Rodrigues = np.array(list(map(float,table.data[column:column+3]))))
    elif inputtype == 'matrix':
      
      o = damask.Orientation(matrix = check_matrix(np.array(list(map(float,table.data[column:column+9]))).reshape(3,3)))
    elif inputtype == 'frame':
      M = np.array(list(map(float,table.data[column[0]:column[0]+3] + \
                                  table.data[column[1]:column[1]+3] + \
                                  table.data[column[2]:column[2]+3]))).reshape(3,3).T
      o = damask.Orientation(matrix = check_matrix(M/np.linalg.norm(M,axis=0)))
    elif inputtype == 'quaternion':
      
      o = damask.Orientation(quaternion = check_quaternion(np.array(list(map(float,table.data[column:column+4])))))

    o.quaternion = r*o.quaternion*R                                                                 # apply additional lab and crystal frame rotations

    for output in options.output:
      if   output == 'quaternion': table.data_append(o.asQuaternion())
      elif output == 'rodrigues':  table.data_append(o.asRodrigues())
      elif output == 'eulers':     table.data_append(o.asEulers(degrees=options.degrees))
      elif output == 'matrix':     table.data_append(o.asMatrix())
      elif output == 'angleaxis':  table.data_append(o.asAngleAxis(degrees=options.degrees,flat=True))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
