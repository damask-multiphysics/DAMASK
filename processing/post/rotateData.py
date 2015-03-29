#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math
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
Rotate vector and/or tensor column data by given angle around given axis.

""", version = scriptID)

parser.add_option('-v','--vector',  dest = 'vector', action = 'extend', metavar = '<string LIST>',
                                    help = 'column heading of vector to rotate')
parser.add_option('-t','--tensor',  dest = 'tensor', action = 'extend', metavar = '<string LIST>',
                                    help = 'column heading of tensor to rotate')
parser.add_option('-r', '--rotation',dest = 'rotation', type = 'float', nargs = 4, metavar = ' '.join(['float']*4),
                                    help = 'angle and axis to rotate data %default')
parser.add_option('-d', '--degrees', dest = 'degrees', action = 'store_true',
                                    help = 'angles are given in degrees [%default]')

parser.set_defaults(vector = [])
parser.set_defaults(tensor = [])
parser.set_defaults(rotation = (0.,1.,1.,1.))                                                       # no rotation about 1,1,1
parser.set_defaults(degrees = False)

(options,filenames) = parser.parse_args()

datainfo = {                                                                                       # list of requested labels per datatype
             'vector':     {'len':3,
                            'label':[]},
             'tensor':     {'len':9,
                            'label':[]},
           }

if options.vector != []: datainfo['vector']['label'] += options.vector
if options.tensor != []: datainfo['tensor']['label'] += options.tensor

toRadians = math.pi/180.0 if options.degrees else 1.0                                               # rescale degrees to radians
r = damask.Quaternion().fromAngleAxis(toRadians*options.rotation[0],options.rotation[1:])
R = r.asMatrix()

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
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out columns to process  ---------------------------------------------------
  active = defaultdict(list)
  column = defaultdict(dict)

  for datatype,info in datainfo.items():
    for label in info['label']:
      key = '1_'+label
      if key in table.labels:
        active[datatype].append(label)
        column[datatype][label] = table.labels.index(key)                                           # remember columns of requested data
      else:
        file['croak'].write('column %s not found...\n'%label)
       
# ------------------------------------------ assemble header ---------------------------------------
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table

    datatype = 'vector'

    for label in active[datatype] if datatype in active else []:                                    # loop over all requested labels
      table.data[column[datatype][label]:column[datatype][label]+datainfo[datatype]['len']] = \
        r * np.array(map(float,
                         table.data[column[datatype][label]:\
                                    column[datatype][label]+datainfo[datatype]['len']]))

    datatype = 'tensor'

    for label in active[datatype] if datatype in active else []:                                    # loop over all requested labels
      A = np.array(map(float,table.data[column[datatype][label]:\
                                        column[datatype][label]+datainfo[datatype]['len']])).\
                  reshape(np.sqrt(datainfo[datatype]['len']),
                          np.sqrt(datainfo[datatype]['len']))
      table.data[column[datatype][label]:\
                 column[datatype][label]+datainfo[datatype]['len']] = \
          np.dot(R,np.dot(A,R.transpose())).reshape(datainfo[datatype]['len'])
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush() 

  table.input_close()                                                                               # close input ASCII table (works for stdin)
  table.output_close()                                                                              # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
