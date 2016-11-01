#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,sys
import math                                                                                         # noqa
import numpy as np
from optparse import OptionParser, OptionGroup
import damask

#"https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions"
def periodicAverage(Points, Box):
  theta = (Points/Box[1]) * (2.0*np.pi)
  xi    = np.cos(theta)
  zeta  = np.sin(theta)
  theta_avg = np.arctan2(-1.0*zeta.mean(), -1.0*xi.mean()) + np.pi
  Pmean = Box[1] * theta_avg/(2.0*np.pi)
  return Pmean

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

group = OptionGroup(parser, "periodic averaging", "")

group.add_option('-p','--periodic',
                  dest = 'periodic',
                  action = 'store_true',
                  help = 'calculate average in periodic space defined by periodic length [%default]')
group.add_option('--boundary',
                 dest = 'boundary', metavar = 'MIN MAX',
                 type = 'float', nargs = 2,
                 help = 'define periodic box end points %default')

parser.add_option_group(group)

parser.set_defaults(function = 'np.average',
                    all      = False,
                    periodic = False,
                    boundary = [0.0, 1.0])

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

  table.data = table.data[np.lexsort([table.data[:,grpColumn]])]                                    # sort data by grpColumn
  
  values,index = np.unique(table.data[:,grpColumn], return_index = True)                            # unique grpColumn values and their positions
  index = np.append(index,rows)                                                                     # add termination position
  grpTable = np.empty((len(values), cols))                                                          # initialize output
  
  for i in range(len(values)):                                                                      # iterate over groups (unique values in grpColumn)
    if options.periodic :
       grpTable[i] = periodicAverage(table.data[index[i]:index[i+1]],options.boundary)              # apply periodicAverage mapping function
    else :
       grpTable[i] = np.apply_along_axis(mapFunction,0,table.data[index[i]:index[i+1]])             # apply mapping function
    if not options.all: grpTable[i,grpColumn] = table.data[index[i],grpColumn]                      # restore grouping column value
  
  table.data = grpTable

# ------------------------------------------ output result -------------------------------  

  table.data_writeArray()
  table.close()                                                                                     # close ASCII table
