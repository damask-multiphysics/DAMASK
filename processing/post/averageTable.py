#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,string
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id: AverageTable.py 3878 2015-02-22 19:26:52Z t.maiti $','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Replace all rows for which the indicator column has identical values by a single row containing their average.
Output table will contain as many rows as there are different (unique) values in the indicator column.

Examples:
For grain averaged values, replace all rows of particular #texture# with a single row containing their average.
""", version = scriptID)

parser.add_option('-l','--label',   dest='key', type="string", metavar='label',
                  help='column label for averaging rows')
(options,filenames) = parser.parse_args()

if options.key == None:
  parser.error('No sorting column specified.')


# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

# ------------------------------------------ loop over input files -----------------------  

for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  table = damask.ASCIItable(file['input'],file['output'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info
  table.info_append(string.replace(scriptID,'\n','\\n') + \
                    '\t' + ' '.join(sys.argv[1:]))

# ------------------------------------------ assemble header -----------------------------  

  table.head_write()

# ------------------------------------------ process data -------------------------------- 

  rows, cols = table.data_readArray()

  table.data = table.data[np.lexsort([table.data[:,table.labels_index(options.key)]])]
  
  values, index = np.unique(table.data[:,table.labels_index(options.key)], return_index=True)
  index = np.append(index,rows)
  avgTable = np.empty((len(values), cols))
  
  for j in xrange(cols) :
    for i in xrange(len(values)) :
      avgTable[i,j] = np.average(table.data[index[i]:index[i+1],j])
  
  table.data = avgTable
  table.data_writeArray()
# ------------------------------------------ output result -------------------------------  

  table.output_flush()                                                      # just in case of buffered ASCII table

  table.input_close()                                                       # close input ASCII table
  if file['name'] != 'STDIN':
    table.output_close()                                                    # close output ASCII table
    os.rename(file['name']+'_tmp',options.key+'_averaged_'+file['name'])                             # overwrite old one with tmp new
