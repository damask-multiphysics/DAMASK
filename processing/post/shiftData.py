#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Uniformly shift column values by given offset.

""", version = scriptID)

parser.add_option('-l','--label',
                  dest = 'label',
                  action = 'extend', metavar = '<string LIST>',
                  help  ='column(s) to shift')
parser.add_option('-o','--offset',
                  dest = 'offset',
                  action = 'extend', metavar='<float LIST>',
                  help = 'offset(s) per column')

parser.set_defaults(label  = [],
                   )

(options,filenames) = parser.parse_args()

if len(options.label) != len(options.offset):
  parser.error('number of column labels and offsets do not match.')

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

  errors  = []
  remarks = []
  columns  = []
  dims     = []
  offsets  = []

  for what,offset in zip(options.label,options.offset):
    col = table.label_index(what)
    if col < 0: remarks.append('column {} not found...'.format(what,type))
    else:
      columns.append(col)
      offsets.append(float(offset))
      dims.append(table.label_dimension(what))

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue
       
# ------------------------------------------ assemble header ---------------------------------------  

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for col,dim,offset in zip(columns,dims,offsets):                                                # loop over items
      table.data[col:col+dim] = offset + np.array(table.data[col:col+dim],'d')
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
