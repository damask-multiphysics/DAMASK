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
Sort rows by given column label(s).

Examples:
With coordinates in columns "x", "y", and "z"; sorting with x slowest and z fastest varying index: --label x,y,z.
""", version = scriptID)


parser.add_option('-l','--label',
                  dest   = 'keys',
                  action = 'extend', metavar = '<string LIST>',
                  help   = 'list of column labels (a,b,c,...)')
parser.add_option('-r','--reverse',
                  dest   = 'reverse',
                  action = 'store_true',
                  help   = 'sort in reverse')

parser.set_defaults(key = [],
                    reverse = False,
                   )

(options,filenames) = parser.parse_args()

if options.keys == None:
  parser.error('No sorting column(s) specified.')

options.keys.reverse()                                                                              # numpy sorts with most significant column as last

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ assemble header ---------------------------------------  

  table.head_read()
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data ---------------------------------------  

  table.data_readArray()
  cols    = []
  remarks = []
  for i,column in enumerate(table.label_index(options.keys)):
    if column < 0:
      remarks.append("label {0} not present.".format(options.keys[i]))
    else:
      cols += [table.data[:,column]]
  if remarks != []: damask.util.croak(remarks)
    
  ind = np.lexsort(cols) if cols != [] else np.arange(table.data.shape[0])
  if options.reverse: ind = ind[::-1]

# ------------------------------------------ output result ---------------------------------------  

  table.data = table.data[ind]
  table.data_writeArray()
  table.close()                                                                                     # close ASCII table
