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
Add column(s) containing eigenvalues and eigenvectors of requested symmetric tensor column(s).

""", version = scriptID)

parser.add_option('-t','--tensor',
                  dest = 'tensor',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'heading of columns containing tensor field values')
parser.add_option('--no-check',
                  dest = 'rh',
                  action = 'store_false',
                  help = 'skip check for right-handed eigenvector basis')

parser.set_defaults(rh = True,
                    )
(options,filenames) = parser.parse_args()

if options.tensor is None:
  parser.error('no data column specified.')

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


# ------------------------------------------ assemble header 1 ------------------------------------

  items = {
            'tensor': {'dim': 9, 'shape': [3,3], 'labels':options.tensor, 'column': []},
          }
  errors  = []
  remarks = []
  
  for type, data in items.items():
    for what in data['labels']:
      dim = table.label_dimension(what)
      if dim != data['dim']: remarks.append('column {} is not a {}...'.format(what,type))
      else:
        items[type]['column'].append(table.label_index(what))
        for order in ['Min','Mid','Max']:
          table.labels_append(['eigval{}({})'.format(order,what)])                                         # extend ASCII header with new labels
        for order in ['Min','Mid','Max']:
          table.labels_append(['{}_eigvec{}({})'.format(i+1,order,what) for i in range(3)])                # extend ASCII header with new labels

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header 2 ------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data -----------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for type, data in items.items():
      for column in data['column']:
        (u,v) = np.linalg.eigh(np.array(list(map(float,table.data[column:column+data['dim']]))).reshape(data['shape']))
        if options.rh and np.dot(np.cross(v[:,0], v[:,1]), v[:,2]) < 0.0 : v[:, 2] *= -1.0          # ensure right-handed eigenvector basis
        table.data_append(list(u))                                                                  # vector of max,mid,min eigval
        table.data_append(list(v.transpose().reshape(data['dim'])))                                 # 3x3=9 combo vector of max,mid,min eigvec coordinates
    outputAlive = table.data_write()                                                                # output processed line in accordance with column labeling

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)
