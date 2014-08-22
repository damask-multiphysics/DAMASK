#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = scriptID.split()[1][:-3]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Permute all values in given column(s).

""", version = scriptID)

parser.add_option('-l','--label',   dest='label', action='extend', metavar='<strig LIST>',
                                    help='heading(s) of column to permute')
parser.add_option('-r', '--rnd',    dest='randomSeed', type='int', metavar='int',
                                    help='seed of random number generator [%default]')

parser.set_defaults(label = [])
parser.set_defaults(randomSeed = None)

(options,filenames) = parser.parse_args()

if len(options.label)== 0:
  parser.error('no data column specified...')

datainfo = {                                                                                        # list of requested labels per datatype
             'scalar':     {'len':1,
                            'label':[]},
           }

if options.label != None: datainfo['scalar']['label'] += options.label
np.random.seed(options.randomSeed)

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

# --------------- figure out columns to process  ---------------------------------------------------
  active = []
  column = {}

  for label in datainfo['scalar']['label']:
    if label in table.labels:
      active.append(label)
      column[label] = table.labels.index(label)                                                     # remember columns of requested data
    else:
      file['croak'].write('column %s not found...\n'%label)
       
# ------------------------------------------ assemble header ---------------------------------------
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  permutation = {}
  table.data_readArray([column[label] for label in active])
  for i,label in enumerate(active):
    mySeed = np.random.randint(9999999)
    np.random.seed(mySeed)
    unique = list(set(table.data[:,i]))
    permutated = np.random.permutation(unique)
    permutation[label] = dict(zip(unique,permutated))

  table.data_rewind()
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for label in active:                                                                            # loop over all requested stiffnesses
      for c in xrange(column[label],
                      column[label]+datainfo['scalar']['len']):
        table.data[c] = permutation[label][float(table.data[c])]                                    # apply permutation
    
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------  
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table
  table.output_close()                                                                              # close output ASCII table
  os.rename(file['name']+'_tmp',file['name'])                                                       # overwrite old one with tmp new
