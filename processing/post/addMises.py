#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,numpy,string,damask
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



def Mises(what,tensor):

  dev = tensor - numpy.trace(tensor)/3.0*numpy.eye(3)
  symdev = 0.5*(dev+dev.T)
  return math.sqrt(numpy.sum(symdev*symdev.T)*
        {
         'stress': 3.0/2.0,
         'strain': 2.0/3.0,
         }[what.lower()])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Add vonMises equivalent values for symmetric part of requested strains and/or stresses.

""" + string.replace(scriptID,'\n','\\n')
)


parser.add_option('-e','--strain',      dest='strain', action='extend', type='string', \
                                        help='heading(s) of columns containing strain tensors')
parser.add_option('-s','--stress',      dest='stress', action='extend', type='string', \
                                        help='heading(s) of columns containing stress tensors')

parser.set_defaults(strain = [])
parser.set_defaults(stress = [])

(options,filenames) = parser.parse_args()

if len(options.strain) + len(options.stress) == 0:
  parser.error('no data column specified...')

datainfo = {                                                               # list of requested labels per datatype
             'strain':     {'len':9,
                            'label':[]},
             'stress':     {'len':9,
                            'label':[]},
           }


if options.strain != None:    datainfo['strain']['label'] += options.strain
if options.stress != None:    datainfo['stress']['label'] += options.stress

# ------------------------------------------ setup file handles ---------------------------------------

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  table = damask.ASCIItable(file['input'],file['output'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info
  table.info_append(string.replace(scriptID,'\n','\\n') + '\t' + ' '.join(sys.argv[1:]))

  active = {}
  column = {}
  head = []

  for datatype,info in datainfo.items():
    for label in info['label']:
      key = {True :'1_%s',
             False:'%s'   }[info['len']>1]%label
      if key not in table.labels:
        sys.stderr.write('column %s not found...\n'%key)
      else:
        if datatype not in active: active[datatype] = []
        if datatype not in column: column[datatype] = {}
        active[datatype].append(label)
        column[datatype][label] = table.labels.index(key)                   # remember columns of requested data
        table.labels_append('Mises(%s)'%label)                                # extend ASCII header with new labels

# ------------------------------------------ assemble header ---------------------------------------  

  table.head_write()

# ------------------------------------------ process data ---------------------------------------  

  while table.data_read():                                                  # read next data line of ASCII table
  
    for datatype,labels in active.items():                                  # loop over vector,tensor
      for label in labels:                                                  # loop over all requested norms
        table.data_append(Mises(datatype,
                                numpy.array(map(float,table.data[column[datatype][label]:
                                                                 column[datatype][label]+datainfo[datatype]['len']]),'d').reshape(3,3)))

    table.data_write()                                                      # output processed line

# ------------------------------------------ output result ---------------------------------------  

  table.output_flush()                                                      # just in case of buffered ASCII table

  file['input'].close()                                                     # close input ASCII table
  if file['name'] != 'STDIN':
    file['output'].close                                                    # close output ASCII table
    os.rename(file['name']+'_tmp',file['name'])                             # overwrite old one with tmp new

