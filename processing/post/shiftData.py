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
Shift values of scalar/special, vector, or tensor columns by given offset.

""", version = scriptID)

parser.add_option('-s','--special',  dest='special', action='extend', metavar='<string LIST>',
                                     help='heading of columns containing field values of special dimension')
parser.add_option('-d','--dimension',dest='N', type='int', metavar='int',
                                     help='dimension of special field values [%default]')
parser.add_option('-v','--vector',   dest='vector', action='extend', metavar='<string LIST>',
                                     help='column heading to shift by vector')
parser.add_option('-t','--tensor',   dest='tensor', action='extend', metavar='<string LIST>',
                                     help='column heading to shift by tensor')
parser.add_option('-o','--offset',   dest='delta', action='extend', metavar='<float LIST>',
                                     help='list of scalar/special, vector, and tensor shifts (in this order!)')

parser.set_defaults(special = [])
parser.set_defaults(vector = [])
parser.set_defaults(tensor = [])
parser.set_defaults(delta  = [])
parser.set_defaults(N = 1)

(options,filenames) = parser.parse_args()

options.delta = np.array(options.delta,'d')
datainfo = {                                                                                        # list of requested labels per datatype
             'special':    {'len':options.N,
                            'label':[]},
             'vector':     {'len':3,
                            'label':[]},
             'tensor':     {'len':9,
                            'label':[]},
           }

length = 0
if options.special != []: datainfo['special']['label'] += options.special; length += len(options.special)
if options.vector  != []: datainfo['vector']['label']  += options.vector;  length += len(options.vector)
if options.tensor  != []: datainfo['tensor']['label']  += options.tensor;  length += len(options.tensor)
if len(options.delta) != length:
  parser.error('length of offset vector does not match column types...')

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
  active = defaultdict(list)
  column = defaultdict(dict)

  for datatype,info in datainfo.items():
    for label in info['label']:
      key = '1_'+label if info['len'] > 1 else label                                                # non-special labels have to start with '1_'
      if key in table.labels:
          active[datatype].append(label)
          column[datatype][label] = table.labels.index(key)                                         # remember columns of requested data
      else:
        file['croak'].write('column %s not found...\n'%label)
       
# ------------------------------------------ assemble header ---------------------------------------
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    i = 0
    for datatype,labels in sorted(active.items(),key=lambda x:datainfo[x[0]]['len']):               # loop over scalar,vector,tensor
      for label in labels:                                                                          # loop over all requested labels
        for j in xrange(datainfo[datatype]['len']):                                                 # loop over entity elements
          table.data[column[datatype][label]+j] = float(table.data[column[datatype][label]+j]) + options.delta[i]
        i += 1  
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table (works for stdin)
  table.output_close()                                                                              # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
