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
Add cumulative (sum of first to current row) values for given label(s).
""", version = scriptID)

parser.add_option('-l','--label',
                  dest='label',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'columns to cumulate')

parser.add_option('-p','--product',
                  dest='product', action = 'store_true',
                  help = 'product of values instead of sum')

(options,filenames) = parser.parse_args()

if options.label is None:
  parser.error('no data column(s) specified.')

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except IOError: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------  

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  errors  = []
  remarks = []
  columns = []
  dims    = []
  
  for what in options.label:
    dim = table.label_dimension(what)
    if dim < 0: remarks.append('column {} not found...'.format(what))
    else:
      dims.append(dim)
      columns.append(table.label_index(what))
      table.labels_append('cum({})'.format(what) if dim == 1 else
                         ['{}_cum({})'.format(i+1,what) for i in range(dim)]  )                     # extend ASCII header with new labels

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header ---------------------------------------  

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data ------------------------------------------ 
  mask = []
  for col,dim in zip(columns,dims): mask += range(col,col+dim)                                      # isolate data columns to cumulate
  cumulated = np.ones(len(mask)) if options.product else np.zeros(len(mask))                        # prepare output field

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    if options.product:
      for i,col in enumerate(mask):
        cumulated[i] *= float(table.data[col])                                                      # cumulate values (multiplication)
    else:
      for i,col in enumerate(mask):
        cumulated[i] += float(table.data[col])                                                      # cumulate values (addition)
    table.data_append(cumulated)

    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
