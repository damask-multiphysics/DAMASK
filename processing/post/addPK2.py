#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,string
import numpy as np
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID = '$Id$'
scriptName = scriptID.split()[1]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing Second Piola--Kirchhoff stress based on given column(s) of
deformation gradient and first Piola--Kirchhoff stress.

""", version = string.replace(scriptID,'\n','\\n')
)

parser.add_option('-f','--defgrad',     dest='defgrad', action='store', type='string', metavar='string',
                                        help='heading of columns containing deformation gradient [%default]')
parser.add_option('-p','--stress',      dest='stress', action='store', type='string', metavar='string',
                                        help='heading of columns containing first Piola--Kirchhoff stress [%default]')
parser.set_defaults(defgrad = 'f')
parser.set_defaults(stress = 'p')

(options,filenames) = parser.parse_args()

datainfo = {                                                                                        # list of requested labels per datatype
             'defgrad':    {'len':9,
                            'label':[]},
             'stress':     {'len':9,
                            'label':[]},
           }

datainfo['defgrad']['label'].append(options.defgrad)
datainfo['stress']['label'].append(options.stress)

# ------------------------------------------ setup file handles ---------------------------------------
files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

# ------------------------------------------ loop over input files ---------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  table = damask.ASCIItable(file['input'],file['output'],False)                                     # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.info_append(string.replace(scriptID,'\n','\\n') + '\t' + ' '.join(sys.argv[1:]))

  active = defaultdict(list)
  column = defaultdict(dict)
  missingColumns = False
  
  for datatype,info in datainfo.items():
    for label in info['label']:
      key = '1_%s'%label
      if key not in table.labels:
        file['croak'].write('column %s not found...\n'%key)
        missingColumns = True                                                                       # break if label not found
      else:
        active[datatype].append(label)
        column[datatype][label] = table.labels.index(key)                                           # remember columns of requested data

  if missingColumns:
    continue

 # ------------------------------------------ assemble header ------------------------------------ 
  table.labels_append(['%i_S'%(i+1) for i in xrange(datainfo['stress']['len'])])                    # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ process data ----------------------------------------  
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    F = np.array(map(float,table.data[column['defgrad'][active['defgrad'][0]]:
                                      column['defgrad'][active['defgrad'][0]]+
                                             datainfo['defgrad']['len']]),'d').reshape(3,3)
    P = np.array(map(float,table.data[column['stress'][active['stress'][0]]:
                                      column['stress'][active['stress'][0]]+
                                             datainfo['stress']['len']]),'d').reshape(3,3)

    table.data_append(list(np.dot(np.linalg.inv(F),P).reshape(9)))                                  # [S] =[P].[F-1]
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result ---------------------------------------  
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  file['input'].close()                                                                             # close input ASCII table (works for stdin)
  file['output'].close()                                                                            # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
