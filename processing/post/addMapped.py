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
Add data in column(s) of second ASCIItable selected from row that is given by the value in a mapping column.

""", version = scriptID)

parser.add_option('-a','--asciitable',  dest='asciitable', metavar='string',
                                        help='mapped ASCIItable')
parser.add_option('-c','--map',         dest='map', metavar='string',
                                        help='heading of column containing row mapping')
parser.add_option('-o','--offset',      dest='offset', type='int', metavar='int',
                                        help='offset between mapped column value and row')
parser.add_option('-v','--vector',      dest='vector', action='extend', metavar='<string LIST>',
                                        help='heading of columns containing vector field values')
parser.add_option('-t','--tensor',      dest='tensor', action='extend', metavar='<string LIST>',
                                        help='heading of columns containing tensor field values')
parser.add_option('-s','--special',     dest='special', action='extend', metavar='<string LIST>',
                                        help='heading of columns containing field values of special dimension')
parser.add_option('-d','--dimension',   dest='N', type='int', metavar='int',
                                        help='dimension of special field values [%default]')
parser.set_defaults(vector = [])
parser.set_defaults(tensor = [])
parser.set_defaults(special = [])
parser.set_defaults(offset = 0)
parser.set_defaults(N = 1)

(options,filenames) = parser.parse_args()

if len(options.vector) + len(options.tensor) + len(options.special) == 0:
  parser.error('no data column specified...')
if options.map == None:
  parser.error('missing mapping column...')

datainfo = {                                                                                        # list of requested labels per datatype
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

# ------------------------------------------ processing mapping ASCIItable -------------------------
if options.asciitable != None and os.path.isfile(options.asciitable):
  mappedTable = damask.ASCIItable(open(options.asciitable),None,False) 
  mappedTable.head_read()                                                                           # read ASCII header info of mapped table

  labels  = []
  indices = []

  for datatype,info in datainfo.items():
    for label in info['label']:
      key = '1_'+label if info['len'] > 1 else label
      if key in mappedTable.labels:
        labels.append(label)                                                                        # extend labels
        indices += range(mappedTable.labels.index(key),
                         mappedTable.labels.index(key)+datainfo[datatype]['len'])
      else:
        sys.stderr.write('column %s not found...\n'%label)
        break

  mappedTable.data_readArray(indices)
  mappedTable.input_close()                                                                        # close mapped input ASCII table

else:
  parser.error('missing mapped ASCIItable...')

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

  if options.map not in table.labels:
    file['croak'].write('column %s not found...\n'%options.map)
    continue

# ------------------------------------------ assemble header --------------------------------------
  for datatype,info in datainfo.items():
    for label in info['label']:
      table.labels_append(label if info['len'] == 1 else \
                          ['%i_%s'%(i+1,label) for i in xrange(info['len'])])                       # extend ASCII header of current table with new labels
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  mappedColumn = table.labels.index(options.map)
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    table.data_append(mappedTable.data[int(table.data[mappedColumn])+options.offset-1])             # add all mapped data types
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table (works for stdin)
  table.output_close()                                                                              # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
