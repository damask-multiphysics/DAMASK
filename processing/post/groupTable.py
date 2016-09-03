#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,sys
import math                                                                                       # noqa
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Apply a user-specified function to condense all rows for which column 'label' has identical values into a single row.
Output table will contain as many rows as there are different (unique) values in the grouping column.

Examples:
For grain averaged values, replace all rows of particular 'texture' with a single row containing their average.
""", version = scriptID)

parser.add_option('-l','--label',
                  dest = 'label',
                  type = 'string', metavar = 'string',
                  help = 'column label for grouping rows')
parser.add_option('-f','--function',
                  dest = 'function',
                  type = 'string', metavar = 'string',
                  help = 'mapping function [%default]')
parser.add_option('-a','--all',
                  dest = 'all',
                  action = 'store_true',
                  help = 'apply mapping function also to grouping column')

parser.set_defaults(function = 'np.average')

(options,filenames) = parser.parse_args()

funcModule,funcName = options.function.split('.')

try:
  mapFunction = getattr(locals().get(funcModule) or 
                        globals().get(funcModule) or
                        __import__(funcModule), 
                        funcName)
except:
  mapFunction = None

if options.label is None:
  parser.error('no grouping column specified.')
if not hasattr(mapFunction,'__call__'):
  parser.error('function "{}" is not callable.'.format(options.function))


# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ sanity checks ---------------------------------------  

  table.head_read()
  if table.label_dimension(options.label) != 1:
    damask.util.croak('column {} is not of scalar dimension.'.format(options.label))
    table.close(dismiss = True)                                                                     # close ASCIItable and remove empty file
    continue
  else:
    grpColumn = table.label_index(options.label)

# ------------------------------------------ assemble info ---------------------------------------  

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data -------------------------------- 

  table.data_readArray()
  rows,cols  = table.data.shape

  table.data = table.data[np.lexsort([table.data[:,grpColumn]])]                                  # sort data by grpColumn
  
  values,index = np.unique(table.data[:,grpColumn], return_index = True)                          # unique grpColumn values and their positions
  index = np.append(index,rows)                                                                   # add termination position
  grpTable = np.empty((len(values), cols))                                                        # initialize output
  
  for i in xrange(len(values)):                                                                   # iterate over groups (unique values in grpColumn)
    grpTable[i] = np.apply_along_axis(mapFunction,0,table.data[index[i]:index[i+1]])              # apply mapping function
    if not options.all: grpTable[i,grpColumn] = table.data[index[i],grpColumn]                    # restore grouping column value
  
  table.data = grpTable

# ------------------------------------------ output result -------------------------------  

  table.data_writeArray()
  table.close()                                                                                     # close ASCII table
