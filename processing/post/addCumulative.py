#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add cumulative (sum of first to current row) values for given label(s).
""", version = scriptID)

parser.add_option('-l','--label',
                  dest='label',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'columns to cumulate')

parser.set_defaults(label = [],
                   )
                    
(options,filenames) = parser.parse_args()

if len(options.label) == 0:
  parser.error('no data column(s) specified.')

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
                         ['{}_cum({})'.format(i+1,what) for i in xrange(dim)]  )                    # extend ASCII header with new labels

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header ---------------------------------------  

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data ------------------------------------------ 

  table.data_readArray()

  mask = []
  for col,dim in zip(columns,dims): mask += range(col,col+dim)                                      # isolate data columns to cumulate

  cumulated = np.zeros((len(table.data),len(mask)))                                                 # prepare output field
  
  for i,values in enumerate(table.data[:,mask]):
    cumulated[i,:] = cumulated[max(0,i-1),:] + values                                               # cumulate values
  
  table.data = np.hstack((table.data,cumulated))

# ------------------------------------------ output result -----------------------------------------

  table.data_writeArray()

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
