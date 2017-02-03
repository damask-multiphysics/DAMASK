#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,sys
import math                                                                                         # noqa
import numpy as np
from optparse import OptionParser, OptionGroup
import damask

def periodicAverage(coords, limits):
  """Centroid in periodic domain, see https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions"""
  theta     = 2.0*np.pi * (coords - limits[0])/(limits[1] - limits[0])
  theta_avg = np.pi + np.arctan2(-np.sin(theta).mean(axis=0), -np.cos(theta).mean(axis=0))
  return limits[0] + theta_avg * (limits[1] - limits[0])/2.0/np.pi

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Apply a user-specified function to condense all rows for which column 'label' has identical values into a single row.
Output table will contain as many rows as there are different (unique) values in the grouping column.
Periodic domain averaging of coordinate values is supported.

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

group.add_option ('-p','--periodic',
                  dest   = 'periodic',
                  action = 'extend', metavar = '<string LIST>',
                  help   = 'coordinate label(s) to average across periodic domain')
group.add_option ('--limits',
                  dest = 'boundary',
                  type = 'float', metavar = 'float float', nargs = 2,
                  help = 'min and max of periodic domain %default')

parser.add_option_group(group)

parser.set_defaults(function = 'np.average',
                    all      = False,
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
  indexrange = table.label_indexrange(options.periodic) if options.periodic is not None else []
  rows,cols  = table.data.shape

  table.data = table.data[np.lexsort([table.data[:,grpColumn]])]                                    # sort data by grpColumn
  
  values,index = np.unique(table.data[:,grpColumn], return_index = True)                            # unique grpColumn values and their positions
  index = np.append(index,rows)                                                                     # add termination position
  grpTable = np.empty((len(values), cols))                                                          # initialize output
  
  for i in range(len(values)):                                                                      # iterate over groups (unique values in grpColumn)
    grpTable[i] = np.apply_along_axis(mapFunction,0,table.data[index[i]:index[i+1]])                # apply (general) mapping function
    grpTable[i,indexrange] = \
      periodicAverage(table.data[index[i]:index[i+1],indexrange],options.boundary)                  # apply periodicAverage mapping function

    if not options.all: grpTable[i,grpColumn] = table.data[index[i],grpColumn]                      # restore grouping column value
  
  table.data = grpTable

# ------------------------------------------ output result -------------------------------  

  table.data_writeArray()
  table.close()                                                                                     # close ASCII table
