#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,string,damask
from optparse import OptionParser, Option

scriptID = '$Id: addNorm.py 3167 2014-06-06 09:43:28Z p.eisenlohr $'
scriptName = scriptID.split()[1]

# -----------------------------
class extendableOption(Option):
# -----------------------------
# used for definition of new option parser action 'extend', which enables to take multiple option arguments
# taken from online tutorial http://docs.python.org/library/optparse.html
  
  ACTIONS = Option.ACTIONS + ("extend",)
  STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
  TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
  ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

  def take_action(self, action, dest, opt, value, values, parser):
    if action == "extend":
      lvalue = value.split(",")
      values.ensure_value(dest, []).extend(lvalue)
    else:
      Option.take_action(self, action, dest, opt, value, values, parser)



# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Add data in column(s) of second ASCIItable selected from row that is given by the value in a mapping column.

""" + string.replace(scriptID,'\n','\\n')
)

parser.add_option('-a','--asciitable',  dest='asciitable', type='string', metavar='FILE',
                                        help='heading of column containing row mapping')
parser.add_option('-c','--map',         dest='map', type='string', metavar='LABEL',
                                        help='heading of column containing row mapping')
parser.add_option('-o','--offset',      dest='offset', type='int', metavar='n',
                                        help='offset between mapped column value and row')
parser.add_option('-v','--vector',      dest='vector', action='extend', type='string', metavar='LABEL',
                                        help='heading of columns containing vector field values')
parser.add_option('-t','--tensor',      dest='tensor', action='extend', type='string', metavar='LABEL',
                                        help='heading of columns containing tensor field values')
parser.add_option('-s','--special',     dest='special', action='extend', type='string', metavar='LABEL',
                                        help='heading of columns containing field values of special dimension')
parser.add_option('-d','--dimension',   dest='N', action='store', type='int', \
                                        help='dimension of special field values [%default]')

parser.set_defaults(vector = [])
parser.set_defaults(tensor = [])
parser.set_defaults(special = [])
parser.set_defaults(offset = 0)
parser.set_defaults(N = 1)

(options,filenames) = parser.parse_args()

if len(options.vector) + len(options.tensor) + len(options.special) == 0:
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

if options.asciitable != None and os.path.isfile(options.asciitable):
  mappedTable = damask.ASCIItable(open(options.asciitable),open(options.asciitable),False) 
  mappedTable.head_read()                                                   # read ASCII header info of mapped table

  labels  = []
  indices = []

  for datatype,info in datainfo.items():
    for label in info['label']:
      foundIt = False
      for key in ['1_'+label,label]:
        if key in mappedTable.labels:
          foundIt = True
          labels.append(label)                                                   # extend labels
          indices += range(mappedTable.labels.index(key),
                           mappedTable.labels.index(key)+datainfo[datatype]['len'])
      if not foundIt:
        file['croak'].write('column %s not found...\n'%label)
        break

  mappedTable.data_readArray(indices)
  mappedTable.input_close()                                                      # close mapped input ASCII table
  mappedTable.output_close()                                                     # close mapped output (same as input) ASCII table

  sys.stderr.write('%s\n'%(mappedTable.data))
else:
  parser.error("Missing mapped ASCIItable")

if options.map == None:
  parser.error("Missing mapping column")


# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  table = damask.ASCIItable(file['input'],file['output'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info
  table.info_append(string.replace(scriptID,'\n','\\n') + '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out columns to process

  if options.map not in table.labels:
    continue

  mappedColumn = table.labels.index(options.map)
  for label in labels:
    table.labels_append(label)                                              # extend ASCII header of current table with new labels

# ------------------------------------------ assemble header ---------------------------------------  

  table.head_write()

# ------------------------------------------ process data ---------------------------------------  

  outputAlive = True

  while outputAlive and table.data_read():                                  # read next data line of ASCII table
    
#    file['croak'].write('%i\n'%(int(table.data[mappedColumn])+options.offset-1))
    table.data_append(mappedTable.data[int(table.data[mappedColumn])+options.offset-1])                # add all mapped data types
    outputAlive = table.data_write()                                        # output processed line

# ------------------------------------------ output result ---------------------------------------  

  outputAlive and table.output_flush()                                      # just in case of buffered ASCII table

  table.input_close()                                                       # close input ASCII table
  if file['name'] != 'STDIN':
    table.output_close()                                                    # close output ASCII table
    os.rename(file['name']+'_tmp',file['name'])                             # overwrite old one with tmp new
