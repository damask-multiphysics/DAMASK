#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = scriptID.split()[1][:-3]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Remove column(s) containing scalar, vectorial, and/or tensorial data.

""", version = scriptID)

parser.add_option('-v','--vector',      dest='vector', action='extend', type='string', metavar='<string LIST>',
                                        help='heading of columns containing 3x1 vector field values')
parser.add_option('-t','--tensor',      dest='tensor', action='extend', type='string', metavar='<string LIST>',
                                        help='heading of columns containing 3x3 tensor field values')
parser.add_option('-s','--special',     dest='special', action='extend', type='string', metavar='<string LIST>',
                                        help='heading of columns containing field values of special dimension')
parser.add_option('-d','--dimension',   dest='N', action='store', type='int', metavar='int',
                                        help='dimension of special field values [%default]')

parser.set_defaults(vector = [])
parser.set_defaults(tensor = [])
parser.set_defaults(special = [])
parser.set_defaults(N = 1)

(options,filenames) = parser.parse_args()

if len(options.vector) + len(options.tensor) + len(options.special)== 0:
  parser.error('no data column specified...')

datainfo = {                                                               # list of requested labels per datatype
             'vector':     {'len':3,
                            'label':[]},
             'tensor':     {'len':9,
                            'label':[]},
             'special':    {'len':options.N,
                            'label':[]},
           }

if options.vector  != None:    datainfo['vector']['label']  += options.vector
if options.tensor  != None:    datainfo['tensor']['label']  += options.tensor
if options.special != None:    datainfo['special']['label'] += options.special

# ------------------------------------------ setup file handles ------------------------------------
files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

#--- loop over input files -------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  table = damask.ASCIItable(file['input'],file['output'],False)                                     # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out columns to delete  ----------------------------------------------------
  columns = []
  for datatype,info in datainfo.items():
    for label in info['label']:
      key = {True :'1_%s',
             False:'%s'   }[info['len']>1]%label
      if key not in table.labels:
        sys.stderr.write('column %s not found...\n'%key)
      else:
        columns.append([table.labels.index(key),info['len']])                                       # remember column and extent of requested data
        if (info['len'] == 1):
            table.labels.remove(label)                                                              # remove single column head
        else:
          for i in xrange(info['len']):
            table.labels.remove('%i_%s'%(i+1,label))                                                # remove multidimensional column head

  columns.sort(key=lambda x:x[0],reverse=True)                                                      # sort from highest column to delete backwards

# ------------------------------------------ assemble header ---------------------------------------
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for col,len in columns:                                                                         # loop over removal candidates
      del table.data[col:col+len]                                                                   # remove each associated entry
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  file['input'].close()                                                                             # close input ASCII table (works for stdin)
  file['output'].close()                                                                            # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
