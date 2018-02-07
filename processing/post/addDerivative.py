#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,sys
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

def derivative(coordinates,what):
  
  result = np.empty_like(what)
  
  # use differentiation by interpolation
  # as described in http://www2.math.umd.edu/~dlevy/classes/amsc466/lecture-notes/differentiation-chap.pdf

  result[1:-1,:] = + what[1:-1,:] * (2.*coordinates[1:-1]-coordinates[:-2]-coordinates[2:]) / \
                     ((coordinates[1:-1]-coordinates[:-2])*(coordinates[1:-1]-coordinates[2:])) \
                   + what[2:,:] * (coordinates[1:-1]-coordinates[:-2]) / \
                     ((coordinates[2:]-coordinates[1:-1])*(coordinates[2:]-coordinates[:-2])) \
                   + what[:-2,:] * (coordinates[1:-1]-coordinates[2:]) / \
                     ((coordinates[:-2]-coordinates[1:-1])*(coordinates[:-2]-coordinates[2:])) \

  result[0,:]    = (what[0,:] - what[1,:]) / \
                   (coordinates[0] - coordinates[1])
  result[-1,:]   = (what[-1,:] - what[-2,:]) / \
                   (coordinates[-1] - coordinates[-2])

  return result
  
# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing numerical derivative of requested column(s) with respect to given coordinates.

""", version = scriptID)

parser.add_option('-c','--coordinates',
                  dest = 'coordinates',
                  type = 'string', metavar='string',
                  help = 'heading of coordinate column')
parser.add_option('-l','--label',
                  dest = 'label',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'heading of column(s) to differentiate')


(options,filenames) = parser.parse_args()

if options.coordinates is None:
  parser.error('no coordinate column specified.')
if options.label is None:
  parser.error('no data column specified.')

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  errors  = []
  remarks = []
  columns = []
  dims    = []
  
  if table.label_dimension(options.coordinates) != 1:
    errors.append('coordinate column {} is not scalar.'.format(options.coordinates))

  for what in options.label:
    dim = table.label_dimension(what)
    if dim < 0: remarks.append('column {} not found...'.format(what))
    else:
      dims.append(dim)
      columns.append(table.label_index(what))
      table.labels_append('d({})/d({})'.format(what,options.coordinates) if dim == 1 else
                         ['{}_d({})/d({})'.format(i+1,what,options.coordinates) for i in range(dim)]  )  # extend ASCII header with new labels

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  table.data_readArray()

  mask = []
  for col,dim in zip(columns,dims): mask += range(col,col+dim)                                      # isolate data columns to differentiate
  
  differentiated = derivative(table.data[:,table.label_index(options.coordinates)].reshape((len(table.data),1)),
                              table.data[:,mask])                                                   # calculate numerical derivative
  
  table.data = np.hstack((table.data,differentiated))

# ------------------------------------------ output result -----------------------------------------

  table.data_writeArray()

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
