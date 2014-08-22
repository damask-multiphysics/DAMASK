#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = scriptID.split()[1][:-3]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing eigenvalues and eigenvectors of requested tensor column(s).

""", version = scriptID)

parser.add_option('-t','--tensor',      dest='tensor', action='extend', type='string', metavar='<string LIST>',
                                        help='heading of columns containing tensor field values')
parser.set_defaults(tensor = [])

(options,filenames) = parser.parse_args()

if len(options.tensor) == 0:
  parser.error('no data column specified...')

datainfo = {                                                                                        # list of requested labels per datatype
             'tensor':     {'len':9,
                            'label':[]},
           }

datainfo['tensor']['label'] += options.tensor

# ------------------------------------------ setup file handles ------------------------------------
files = []
for name in filenames:
  if os.path.exists(name):
    files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

#--- loop over input files -------------------------------------------------------------------------
for file in files:
  file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],True)                                      # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

  active = []
  column = defaultdict(dict)
  
  for label in datainfo['tensor']['label']:
    key = '1_%s'%label
    if key not in table.labels:
      file['croak'].write('column %s not found...\n'%key)
    else:
      active.append(label)
      column[label] = table.labels.index(key)                                                       # remember columns of requested data

# ------------------------------------------ assemble header ---------------------------------------
  for labels in active: 
    table.labels_append(['%i_eigval(%s)'%(i+1,label) for i in xrange(3)])                           # extend ASCII header with new labels
    table.labels_append(['%i_eigvec(%s)'%(i+1,label) for i in xrange(9)])                           # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for labels in active:                                                                           # loop over requested data
      tensor = np.array(map(float,table.data[column[label]:column[label]+datainfo['tensor']['len']])).\
                                                            reshape((datainfo['tensor']['len']//3,3))
      (u,v) = np.linalg.eig(tensor)
      table.data_append(list(u))
      table.data_append(list(v.transpose().reshape(datainfo['tensor']['len'])))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table (works for stdin)
  table.output_close()                                                                              # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
