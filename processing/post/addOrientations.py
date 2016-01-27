#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math
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
Orientation is given by quaternion, Euler angles, rotation matrix, or crystal frame coordinates (i.e. component vectors of rotation matrix).

""", version = scriptID)

outputChoices = ['quaternion','rodrigues','eulers']
parser.add_option('-o', '--output',
                  dest = 'output',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'output orientation formats {%s}'%(','.join(outputChoices)))
parser.add_option('-r', '--rotation',
                  dest='rotation',
                  type = 'float', nargs = 4, metavar = ' '.join(['float']*4),
                  help = 'angle and axis to (pre)rotate orientation')

parser.add_option('-s', '--symmetry',
                  dest = 'symmetry',
                  type = 'choice', choices = damask.Symmetry.lattices[1:], metavar='string',
                  help = 'crystal symmetry [%default] {{{}}} '.format(', '.join(damask.Symmetry.lattices[1:])))
parser.add_option('-e', '--eulers',
                  dest = 'eulers',
                  type = 'string', metavar = 'string',
                  help = 'Euler angles label')
parser.add_option('-d', '--degrees',
                  dest = 'degrees',
                  action = 'store_true',
                  help = 'Euler angles are given in degrees [%default]')
parser.add_option('-m', '--matrix',
                  dest = 'matrix',
                  type = 'string', metavar = 'string',
                  help = 'orientation matrix label')
parser.add_option('-a',
                  dest = 'a',
                  type = 'string', metavar = 'string',
                  help = 'crystal frame a vector label')
parser.add_option('-b',
                  dest = 'b',
                  type = 'string', metavar = 'string',
                  help = 'crystal frame b vector label')
parser.add_option('-c',
                  dest = 'c',
                  type = 'string', metavar = 'string',
                  help = 'crystal frame c vector label')
parser.add_option('-q', '--quaternion',
                  dest = 'quaternion',
                  type = 'string', metavar = 'string',
                  help = 'quaternion label')

parser.set_defaults(output = [],
                    symmetry = damask.Symmetry.lattices[-1],
                    rotation = (0.,1.,1.,1.),                                                       # no rotation about 1,1,1
                    degrees = False,
                   )

(options, filenames) = parser.parse_args()

options.output = map(lambda x: x.lower(), options.output)
if options.output == [] or (not set(options.output).issubset(set(outputChoices))):
  parser.error('output must be chosen from {}.'.format(', '.join(outputChoices)))

input = [options.eulers     != None,
         options.a          != None and \
         options.b          != None and \
         options.c          != None,
         options.matrix     != None,
         options.quaternion != None,
        ]

if np.sum(input) != 1: parser.error('needs exactly one input format.')

(label,dim,inputtype) = [(options.eulers,3,'eulers'),
                         ([options.a,options.b,options.c],[3,3,3],'frame'),
                         (options.matrix,9,'matrix'),
                         (options.quaternion,4,'quaternion'),
                        ][np.where(input)[0][0]]                                                    # select input label that was requested
toRadians = math.pi/180.0 if options.degrees else 1.0                                               # rescale degrees to radians
r = damask.Quaternion().fromAngleAxis(toRadians*options.rotation[0],options.rotation[1:])           # pre-rotation

# --- loop over input files ------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks -----------------------------------------

  errors  = []
  remarks = []
  
  if not np.all(table.label_dimension(label) == dim):  errors.append('input {} has wrong dimension {}.'.format(label,dim))
  else:  column = table.label_index(label)

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  for output in options.output:
    if output == 'quaternion': table.labels_append(['{}_quat({})'.format(     i+1,options.symmetry) for i in xrange(4)])
    if output == 'rodrigues':  table.labels_append(['{}_rodrigues({})'.format(i+1,options.symmetry) for i in xrange(3)])
    if output == 'eulers':     table.labels_append(['{}_eulers({})'.format(   i+1,options.symmetry) for i in xrange(3)])
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    if inputtype == 'eulers':
      o = damask.Orientation(Eulers   = np.array(map(float,table.data[column:column+3]))*toRadians,
                             symmetry = options.symmetry).reduced()
    elif inputtype == 'matrix':
      o = damask.Orientation(matrix   = np.array(map(float,table.data[column:column+9])).reshape(3,3).transpose(),
                             symmetry = options.symmetry).reduced()
    elif inputtype == 'frame':
      o = damask.Orientation(matrix = np.array(map(float,table.data[column[0]:column[0]+3] + \
                                                         table.data[column[1]:column[1]+3] + \
                                                         table.data[column[2]:column[2]+3])).reshape(3,3),
                             symmetry = options.symmetry).reduced()
    elif inputtype == 'quaternion':
      o = damask.Orientation(quaternion = np.array(map(float,table.data[column:column+4])),
                             symmetry   = options.symmetry).reduced()

    o.quaternion = r*o.quaternion

    for output in options.output:
      if output == 'quaternion':  table.data_append(o.asQuaternion())
      if output == 'rodrigues':   table.data_append(o.asRodrigues())
      if output == 'eulers':      table.data_append(o.asEulers('Bunge', degrees=options.degrees))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
