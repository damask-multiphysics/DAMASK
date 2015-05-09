#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math,string
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# definition of element-wise p-norms for matrices
def normAbs(object):        # p = 1
  return sum(map(abs, object))

def normFrobenius(object):  # p = 2
  return math.sqrt(sum([x*x for x in object]))

def normMax(object):        # p = infinity
  return max(map(abs, object))

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing norm of requested column(s) being either vectors or tensors.

""", version = scriptID)

normChoices = ['abs','frobenius','max']
parser.add_option('-n','--norm',        dest='norm', type='choice', choices=normChoices, metavar='string',
                  help='type of element-wise p-norm [frobenius] {%s}'%(','.join(map(str,normChoices))))
parser.add_option('-v','--vector',      dest='vector', action='extend', metavar='<string LIST>',
                  help='heading of columns containing vector field values')
parser.add_option('-t','--tensor',      dest='tensor', action='extend', metavar='<string LIST>',
                  help='heading of columns containing tensor field values')
parser.add_option('-s','--special',     dest='special', action='extend', metavar='<string LIST>',
                  help='heading of columns containing field values of special dimension')
parser.add_option('-d','--dimension',   dest='N', type='int', metavar='int',
                  help='dimension of special field values [%default]')
parser.set_defaults(norm = 'frobenius')
parser.set_defaults(N = 12)

(options,filenames) = parser.parse_args()

if (not None) in [options.vector,options.tensor,options.special]:
  parser.error('no data column specified...')

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

  table = damask.ASCIItable(file['input'],file['output'],False)                                    # make unbuffered ASCII_table
  table.head_read()                                                                                # read ASCII header info
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

  active = defaultdict(list)
  column = defaultdict(dict)

  for datatype,info in datainfo.items():
    for label in info['label']:
      key = '1_'+label if info['len'] > 1 else label                                               # columns of non-scalars need to start with '1_'
      if key not in table.labels:
        file['croak'].write('column %s not found...\n'%key)
      else:
        active[datatype].append(label)
        column[datatype][label] = table.labels.index(key)                                           # remember columns of requested data

# ------------------------------------------ assemble header ---------------------------------------
  for datatype,labels in active.items():                                                            # loop over vector,tensor
    for label in labels:                                                                            # loop over all requested determinants
      table.labels_append('norm%s(%s)'%(options.norm.capitalize(),label))                           # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for datatype,labels in active.items():                                                          # loop over vector,tensor
      for label in labels:                                                                          # loop over all requested norms
        eval("table.data_append(norm%s(map(float,table.data[column[datatype][label]:"\
                     "column[datatype][label]+datainfo[datatype]['len']])))"%options.norm.capitalize())
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table (works for stdin)
  table.output_close()                                                                              # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
