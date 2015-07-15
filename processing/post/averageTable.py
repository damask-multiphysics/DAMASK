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

parser.add_option('-l','--label',   dest='label', type="string", metavar='string',
                  help='column label for grouping rows')
(options,filenames) = parser.parse_args()

if options.label == None:
  parser.error('No sorting column specified.')


# --- loop over input files -------------------------------------------------------------------------

if filenames == []:
  filenames = ['STDIN']

for name in filenames:
  if name == 'STDIN':
    file = {'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m\n')
  else:
    if not os.path.exists(name): continue
    file = {'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],buffered=False)                            # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info

  if table.label_dimension(options.label) != 1:
    file['croak'].write('column {0} is not of scalar dimension...\n'.format(options.label))
    table.close(dismiss = True)                                                                     # close ASCIItable and remove empty file
    continue


# ------------------------------------------ assemble header -----------------------------  

  table.info_append(string.replace(scriptID,'\n','\\n') + \
                    '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data -------------------------------- 

  table.data_readArray()
  rows,cols  = table.data.shape

  table.data = table.data[np.lexsort([table.data[:,table.label_index(options.label)]])]
  
  values,index = np.unique(table.data[:,table.label_index(options.label)], return_index=True)
  index = np.append(index,rows)
  avgTable = np.empty((len(values), cols))
  
  for j in xrange(cols) :
    for i in xrange(len(values)) :
      avgTable[i,j] = np.average(table.data[index[i]:index[i+1],j])
  
  table.data = avgTable

# ------------------------------------------ output result -------------------------------  

  table.data_writeArray()
  table.output_flush()                                                                               # just in case of buffered ASCII table

  table.close()                                                                                      # close ASCII table
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',options.key+'_averaged_'+file['name'])                             # overwrite old one with tmp new
