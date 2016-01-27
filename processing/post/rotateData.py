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
Rotate vector and/or tensor column data by given angle around given axis.

""", version = scriptID)

parser.add_option('-v','--vector',
                  dest = 'vector',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'column heading of vector(s) to rotate')
parser.add_option('-t','--tensor',
                  dest = 'tensor',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'column heading of tensor(s) to rotate')
parser.add_option('-r', '--rotation',
                  dest = 'rotation',
                  type = 'float', nargs = 4, metavar = ' '.join(['float']*4),
                  help = 'angle and axis to rotate data [%default]')
parser.add_option('-d', '--degrees',
                  dest = 'degrees',
                  action = 'store_true',
                  help = 'angles are given in degrees [%default]')

parser.set_defaults(rotation = (0.,1.,1.,1.),                                                       # no rotation about 1,1,1
                    degrees = False,
                   )
                    
(options,filenames) = parser.parse_args()

if options.vector == None and options.tensor == None:
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

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  items = {
            'tensor': {'dim': 9, 'shape': [3,3], 'labels':options.tensor, 'active':[], 'column': []},
            'vector': {'dim': 3, 'shape': [3],   'labels':options.vector, 'active':[], 'column': []},
          }
  errors  = []
  remarks = []
  column = {}
  
  for type, data in items.iteritems():
    for what in data['labels']:
      dim = table.label_dimension(what)
      if dim != data['dim']: remarks.append('column {} is not a {}.'.format(what,type))
      else:
        items[type]['active'].append(what)
        items[type]['column'].append(table.label_index(what))

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

    datatype = 'vector'

    for column in items[datatype]['column']:                                                        # loop over all requested labels
      table.data[column:column+items[datatype]['dim']] = \
        q * np.array(map(float,table.data[column:column+items[datatype]['dim']]))

    datatype = 'tensor'

    for column in items[datatype]['column']:                                                        # loop over all requested labels
      table.data[column:column+items[datatype]['dim']] = \
        np.dot(R,np.dot(np.array(map(float,table.data[column:column+items[datatype]['dim']])).\
                        reshape(items[datatype]['shape']),R.transpose())).\
                  reshape(items[datatype]['dim'])

    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
