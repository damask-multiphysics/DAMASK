#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing Second Piola--Kirchhoff stress based on given column(s) of deformation gradient and first Piola--Kirchhoff stress.

""", version = scriptID)

parser.add_option('-f','--defgrad',     dest='defgrad', metavar='string',
                  help='heading of columns containing deformation gradient [%default]')
parser.add_option('-p','--stress',      dest='stress', metavar='string',
                  help='heading of columns containing first Piola--Kirchhoff stress [%default]')
parser.set_defaults(defgrad = 'f')
parser.set_defaults(stress  = 'p')

(options,filenames) = parser.parse_args()

# ------------------------------------------ setup file handles ------------------------------------
files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

# ------------------------------------------ loop over input files ---------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  table = damask.ASCIItable(file['input'],file['output'],False)                                     # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out columns to process  ---------------------------------------------------
  missingColumns = False
  column={ 'defgrad': table.labels.index('1_'+options.defgrad), 
           'stress':  table.labels.index('1_'+options.stress)}
  for key in column:
    if column[key]<1:
      file['croak'].write('column %s not found...\n'%key)
      missingColumns=True
  if missingColumns: continue

# ------------------------------------------ assemble header --------------------------------------
  table.labels_append(['%i_S'%(i+1) for i in xrange(9)])                                       # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    F = np.array(map(float,table.data[column['defgrad']:column['defgrad']+9]),'d').reshape(3,3)
    P = np.array(map(float,table.data[column['stress'] :column['stress']+9]),'d').reshape(3,3)
    table.data_append(list(np.dot(np.linalg.inv(F),P).reshape(9)))                                  # [S] =[P].[F-1]
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table (works for stdin)
  table.output_close()                                                                              # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
