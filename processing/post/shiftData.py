#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,numpy,string
import damask
from collections import defaultdict
from optparse import OptionParser, Option

scriptID = '$Id$'
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
Shift values of scalar, vector, or tensor columns by given offset.

""" + string.replace(scriptID,'\n','\\n')
)

parser.add_option('-s','--scalar',  dest='scalar', action='extend', type='string',
                                    help='column heading to shift by scalar',
                                    metavar='<label(s)>')
parser.add_option('-v','--vector',  dest='vector', action='extend', type='string',
                                    help='column heading to shift by vector',
                                    metavar='<label(s)>')
parser.add_option('-t','--tensor',  dest='tensor', action='extend', type='string',
                                    help='column heading to shift by tensor',
                                    metavar='<label(s)>')
parser.add_option('-d','--delta',   dest='delta', action='extend', type='string',
                                    help='list of scalar, vector, and tensor shifts (in this order!)',
                                    metavar='<float(s)>')

parser.set_defaults(scalar = [])
parser.set_defaults(vector = [])
parser.set_defaults(tensor = [])
parser.set_defaults(delta  = [])

(options,filenames) = parser.parse_args()

options.delta = numpy.array(map(float,options.delta))
datainfo = {                                                               # list of requested labels per datatype
             'scalar':     {'len':1,
                            'label':[]},
             'vector':     {'len':3,
                            'label':[]},
             'tensor':     {'len':9,
                            'label':[]},
           }

length = 0
if options.scalar != []: datainfo['scalar']['label'] += options.scalar; length += len(options.scalar)*datainfo['scalar']['len']
if options.vector != []: datainfo['vector']['label'] += options.vector; length += len(options.vector)*datainfo['vector']['len']
if options.tensor != []: datainfo['tensor']['label'] += options.tensor; length += len(options.tensor)*datainfo['tensor']['len']
if len(options.delta) != length:
  parser.error('Length of offset vector does not match column types.')

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
  table.info_append(string.replace(scriptName,'\n','\\n') + \
                    '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out columns to process ---------------------------------------
  active = defaultdict(list)
  column = defaultdict(dict)

  for datatype,info in datainfo.items():
    for label in info['label']:
      foundIt = False
      for key in ['1_'+label,label]:
        if key in table.labels:
          foundIt = True
          active[datatype].append(label)
          column[datatype][label] = table.labels.index(key)                   # remember columns of requested data
      if not foundIt:
        file['croak'].write('column %s not found...\n'%label)
       
# ------------------------------------------ assemble header ---------------------------------------

  table.head_write()

# ------------------------------------------ process data ---------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():                                  # read next data line of ASCII table

    i = 0
    for datatype,labels in sorted(active.items(),key=lambda x:datainfo[x[0]]['len']):  # loop over scalar,vector,tensor
#      file['croak'].write('%s\n'%datatype)
      for label in labels:                                                  # loop over all requested labels
#        file['croak'].write('%s\n'%label)
        for j in xrange(datainfo[datatype]['len']):                         # loop over entity elements
#          file['croak'].write('%s\n'%table.data[column[datatype][label]+j])
          table.data[column[datatype][label]+j] = float(table.data[column[datatype][label]+j]) + options.delta[i]
          i += 1
    
    outputAlive = table.data_write()                                        # output processed line

# ------------------------------------------ output result ---------------------------------------  

  outputAlive and table.output_flush()                                      # just in case of buffered ASCII table

  if file['name'] != 'STDIN':
    file['input'].close()                                                   # close input ASCII table
    file['output'].close()                                                  # close output ASCII table
    os.rename(file['name']+'_tmp',file['name'])                             # overwrite old one with tmp new
