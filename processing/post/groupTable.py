#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser, OptionGroup
import math                                                                                         # noqa

import numpy as np

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

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Apply a user-specified function to condense into a single row all those rows for which columns 'label' have identical values.
Output table will contain as many rows as there are different (unique) values in the grouping column(s).
Periodic domain averaging of coordinate values is supported.

Examples:
For grain averaged values, replace all rows of particular 'texture' with a single row containing their average.
{name} --label texture --function np.average data.txt
""".format(name = scriptName), version = scriptID)

parser.add_option('-l','--label',
                  dest = 'label',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'column label(s) for grouping rows')
parser.add_option('-f','--function',
                  dest = 'function',
                  type = 'string', metavar = 'string',
                  help = 'mapping function [%default]')
parser.add_option('-a','--all',
                  dest = 'all',
                  action = 'store_true',
                  help = 'apply mapping function also to grouping column(s)')
                  
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
                    label    = [],
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

if options.label is []:
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

  remarks = []
  errors = []

  table.head_read()
  grpColumns = table.label_index(options.label)[::-1]
  grpColumns = grpColumns[np.where(grpColumns>=0)]

  if len(grpColumns) == 0:  errors.append('no valid grouping column present.')

  if remarks != []: damask.util.croak(remarks)
  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss=True)
    continue

# ------------------------------------------ assemble info ---------------------------------------  

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data -------------------------------- 

  table.data_readArray()
  indexrange = table.label_indexrange(options.periodic) if options.periodic is not None else []
  rows,cols  = table.data.shape

  table.data = table.data[np.lexsort(table.data[:,grpColumns].T)]                                   # sort data by grpColumn(s)
  values,index = np.unique(table.data[:,grpColumns], axis=0, return_index=True)                     # unique grpColumn values and their positions
  index = sorted(np.append(index,rows))                                                             # add termination position
  grpTable = np.empty((len(values), cols))                                                          # initialize output
  
  for i in range(len(values)):                                                                      # iterate over groups (unique values in grpColumn)
    grpTable[i] = np.apply_along_axis(mapFunction,0,table.data[index[i]:index[i+1]])                # apply (general) mapping function
    grpTable[i,indexrange] = \
      periodicAverage(table.data[index[i]:index[i+1],indexrange],options.boundary)                  # apply periodicAverage mapping function

    if not options.all: grpTable[i,grpColumns] = table.data[index[i],grpColumns]                    # restore grouping column value
  
  table.data = grpTable

# ------------------------------------------ output result -------------------------------  

  table.data_writeArray()
  table.close()                                                                                     # close ASCII table
