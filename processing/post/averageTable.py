#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,string
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Replace all rows for which column 'label' has identical values by a single row containing their average.
Output table will contain as many rows as there are different (unique) values in the grouping column.

Examples:
For grain averaged values, replace all rows of particular 'texture' with a single row containing their average.
""", version = scriptID)

parser.add_option('-l','--label',
                  dest = 'label',
                  type = 'string', metavar = 'string',
                  help = 'column label for grouping rows')

(options,filenames) = parser.parse_args()

if options.label == None:
  parser.error('no grouping column specified.')


# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              outname = options.label+'_averaged_'+name if name else name,
                              buffered = False)
  except: continue
  table.report_name(scriptName,name)

# ------------------------------------------ sanity checks ---------------------------------------  

  table.head_read()
  if table.label_dimension(options.label) != 1:
    table.croak('column {} is not of scalar dimension.'.format(options.label))
    table.close(dismiss = True)                                                                     # close ASCIItable and remove empty file
    continue

# ------------------------------------------ assemble info ---------------------------------------  

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data -------------------------------- 

  table.data_readArray()
  rows,cols  = table.data.shape

  table.data = table.data[np.lexsort([table.data[:,table.label_index(options.label)]])]
  
  values,index = np.unique(table.data[:,table.label_index(options.label)], return_index = True)
  index = np.append(index,rows)
  avgTable = np.empty((len(values), cols))
  
  for j in xrange(cols) :
    for i in xrange(len(values)) :
      avgTable[i,j] = np.average(table.data[index[i]:index[i+1],j])
  
  table.data = avgTable

# ------------------------------------------ output result -------------------------------  

  table.data_writeArray()
  table.close()                                                                                     # close ASCII table
