#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math,string
import numpy as np
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

def Mises(what,tensor):

  dev = tensor - np.trace(tensor)/3.0*np.eye(3)
  symdev = 0.5*(dev+dev.T)
  return math.sqrt(np.sum(symdev*symdev.T)*
        {
         'stress': 3.0/2.0,
         'strain': 2.0/3.0,
         }[what.lower()])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add vonMises equivalent values for symmetric part of requested strains and/or stresses.

""", version = scriptID)

parser.add_option('-e','--strain', dest='strain', action='extend', metavar='<string LIST>',
                  help='heading(s) of columns containing strain tensors')
parser.add_option('-s','--stress', dest='stress', action='extend', metavar='<string LIST>',
                  help='heading(s) of columns containing stress tensors')

(options,filenames) = parser.parse_args()

if (not None) in [options.strain,options.stress]:
  parser.error('no data column specified...')

datainfo = {                                                                                        # list of requested labels per datatype
             'strain':     {'len':9,
                            'label':[]},
             'stress':     {'len':9,
                            'label':[]},
           }

if options.strain != None:    datainfo['strain']['label'] += options.strain
if options.stress != None:    datainfo['stress']['label'] += options.stress

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

  active = defaultdict(list)
  column = defaultdict(dict)

  for datatype,info in datainfo.items():
    for label in info['label']:
      key = '1_%s'%label
      if key not in table.labels:
        file['croak'].write('column %s not found...\n'%key)
      else:
        active[datatype].append(label)
        column[datatype][label] = table.labels.index(key)                                           # remember columns of requested data

# ------------------------------------------ assemble header ---------------------------------------
  for datatype,labels in active.items():                                                            # loop over vector,tensor
    for label in labels:                                                                            # loop over all requested determinants
      table.labels_append('Mises(%s)'%label)                                                        # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for datatype,labels in active.items():                                                          # loop over vector,tensor
      for label in labels:                                                                          # loop over all requested norms
        table.data_append(Mises(datatype,
                    np.array(map(float,table.data[column[datatype][label]:
                                                  column[datatype][label]+datainfo[datatype]['len']]),'d').reshape(3,3)))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table (works for stdin)
  table.output_close()                                                                              # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
