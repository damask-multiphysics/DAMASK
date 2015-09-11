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
Permute all values in given column(s).

""", version = scriptID)

parser.add_option('-l','--label',
                  dest = 'label',
                  action = 'extend', metavar = '<string LIST>',
                  help  ='column(s) to permute')
parser.add_option('-r', '--rnd',
                  dest = 'randomSeed',
                  type = 'int', metavar = 'int',
                  help = 'seed of random number generator [%default]')

parser.set_defaults(label = [],
                    randomSeed = None,
                   )

(options,filenames) = parser.parse_args()

if len(options.label) == 0:
  parser.error('no labels specified.')

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except: continue
  table.report_name(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ process labels ---------------------------------------  

  errors  = []
  remarks = []
  columns = []
  dims    = []

  indices    = table.label_index    (options.label)
  dimensions = table.label_dimension(options.label)
  for i,index in enumerate(indices):
    if index == -1: remarks.append('label {} not present...'.format(options.label[i]))
    else:
      columns.append(index)
      dims.append(dimensions[i])

  if remarks != []: table.croak(remarks)
  if errors  != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue
       
# ------------------------------------------ assemble header ---------------------------------------

  randomSeed = int(os.urandom(4).encode('hex'), 16) if options.randomSeed == None else options.randomSeed         # random seed per file
  np.random.seed(randomSeed)

  table.info_append([scriptID + '\t' + ' '.join(sys.argv[1:]),
                     'random seed {}'.format(randomSeed),
                    ])
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  table.data_readArray()                                                                            # read all data at once
  for col,dim in zip(columns,dims):
    table.data[:,col:col+dim] = np.random.permutation(table.data[:,col:col+dim])

# ------------------------------------------ output result -----------------------------------------  

  table.data_writeArray()

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
