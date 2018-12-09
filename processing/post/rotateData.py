#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Rotate vector and/or tensor column data by given angle around given axis.

""", version = scriptID)

parser.add_option('-d', '--data',
                  dest = 'data',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'vector/tensor value(s) label(s)')
parser.add_option('-r', '--rotation',
                  dest = 'rotation',
                  type = 'float', nargs = 4, metavar = ' '.join(['float']*4),
                  help = 'angle and axis to rotate data [%default]')
parser.add_option('--degrees',
                  dest = 'degrees',
                  action = 'store_true',
                  help = 'angles are given in degrees [%default]')

parser.set_defaults(rotation = (0.,1.,1.,1.),                                                       # no rotation about 1,1,1
                    degrees = False,
                   )
                    
(options,filenames) = parser.parse_args()

if options.data is None:
  parser.error('no data column specified.')

toRadians = math.pi/180.0 if options.degrees else 1.0                                               # rescale degrees to radians
q = damask.Quaternion().fromAngleAxis(toRadians*options.rotation[0],options.rotation[1:])
R = q.asMatrix()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()

  errors  = []
  remarks = []
  active  = {'vector':[],'tensor':[]}

  for i,dim in enumerate(table.label_dimension(options.data)):
    label = options.data[i]
    if dim == -1:
      remarks.append('{} "{}" not found...'.format(datatype,me))
    elif dim ==  3:
      remarks.append('adding vector "{}"...'.format(label))
      active['vector'].append(label)
    elif dim ==  9:
      remarks.append('adding tensor "{}"...'.format(label))
      active['tensor'].append(label)

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for v in active['vector']:
      column = table.label_index(v)
      table.data[column:column+3] = q * np.array(list(map(float,table.data[column:column+3])))
    for t in active['tensor']:
      column = table.label_index(t)
      table.data[column:column+9] = \
        np.dot(R,np.dot(np.array(list(map(float,table.data[column:column+9]))).reshape((3,3)),
                        R.transpose())).reshape((9))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
