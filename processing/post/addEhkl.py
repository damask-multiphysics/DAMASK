#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = scriptID.split()[1][:-3]

def normalize(vec):
    return vec/np.sqrt(np.inner(vec,vec))

def E_hkl(stiffness,vec):   # stiffness = (c11,c12,c44)
    v = normalize(vec)
    S11 = (stiffness[0]+stiffness[1])/(stiffness[0]*stiffness[0]+stiffness[0]*stiffness[1]-2.0*stiffness[1]*stiffness[1])
    S12 = (            -stiffness[1])/(stiffness[0]*stiffness[0]+stiffness[0]*stiffness[1]-2.0*stiffness[1]*stiffness[1])
    S44 = 1.0/stiffness[2]

    invE = S11-(S11-S12-0.5*S44)* (1.0 - \
                 (v[0]**4+v[1]**4+v[2]**4) \
            /#------------------------------------
                 np.inner(v,v)**2 \
                )

    return 1.0/invE

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing directional stiffness based on given cubic stiffness values C11, C12, and C44 in consecutive columns.

""", version = scriptID)

parser.add_option('-c','--stiffness',   dest='vector', action='extend', type='string', metavar='<string LIST>',
                                        help='heading of column containing C11 (followed by C12, C44) field values')
parser.add_option('-d','--direction', \
                       '--hkl',         dest='hkl', action='store', type='int', nargs=3, metavar='int int int',
                                        help='direction of elastic modulus %default')
parser.set_defaults(vector = [])
parser.set_defaults(hkl = [1,1,1])

(options,filenames) = parser.parse_args()

if len(options.vector)== 0:
  parser.error('no data column specified...')

datainfo = {                                                                                        # list of requested labels per datatype
             'vector':     {'len':3,
                            'label':[]},
           }

datainfo['vector']['label']  += options.vector

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

  active = []
  column = defaultdict(dict)

  for label in datainfo['vector']['label']:
    key = '1_%s'%label
    if key not in table.labels:
      file['croak'].write('column %s not found...\n'%key)
    else:
      active.append(label)
      column[label] = table.labels.index(key)                                                       # remember columns of requested data

# ------------------------------------------ assemble header ---------------------------------------
  for label in active:
    table.labels_append('E%i%i%i(%s)'%(options.hkl[0],
                                       options.hkl[1],
                                       options.hkl[2],label))                                       # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for label in active:
      table.data_append(E_hkl(map(float,table.data[column[label]:\
                                                   column[label]+datainfo['vector']['len']]),options.hkl))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table (works for stdin)
  table.output_close()                                                                              # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
