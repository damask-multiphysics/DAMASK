#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys
import numpy as np
from optparse import OptionParser
import damask

scriptID = '$Id$'
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Uniformly scale column values by given factor.

""", version = scriptID)

parser.add_option('-l','--label',
                  dest = 'label',
                  action = 'extend', metavar = '<string LIST>',
                  help  ='column(s) to scale')
parser.add_option('-f','--factor',
                  dest = 'factor',
                  action = 'extend', metavar='<float LIST>',
                  help = 'factor(s) per column')

parser.set_defaults(label  = [],
                   )

(options,filenames) = parser.parse_args()

if len(options.label) != len(options.factor):
  parser.error('number of column labels and factors do not match.')

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
  factors  = []

  for what,factor in zip(options.label,options.factor):
    col = table.label_index(what)
    if col < 0: remarks.append('column {} not found...'.format(what,type))
    else:
      columns.append(col)
      factors.append(float(factor))
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
    for col,dim,factor in zip(columns,dims,factors):                                                # loop over items
      table.data[col:col+dim] = factor * np.array(table.data[col:col+dim],'d')
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
