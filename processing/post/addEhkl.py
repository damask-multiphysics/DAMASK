#!/usr/bin/env python

import os,re,sys,math,numpy,string,damask
from optparse import OptionParser, Option

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



def normalize(vec):
    return vec/numpy.sqrt(numpy.inner(vec,vec))

def E_hkl(stiffness,vec):   # stiffness = (c11,c12,c44)
    v = normalize(vec)
    S11 = (stiffness[0]+stiffness[1])/(stiffness[0]*stiffness[0]+stiffness[0]*stiffness[1]-2.0*stiffness[1]*stiffness[1])
    S12 = (            -stiffness[1])/(stiffness[0]*stiffness[0]+stiffness[0]*stiffness[1]-2.0*stiffness[1]*stiffness[1])
    S44 = 1.0/stiffness[2]

    invE = S11-(S11-S12-0.5*S44)* (1.0 - \
                 (v[0]**4+v[1]**4+v[2]**4) \
            /#------------------------------------
                 numpy.inner(v,v)**2 \
                )

    return 1.0/invE

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing directional stiffness
based on given cubic stiffness values C11, C12, and C44 in consecutive columns.

""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-c','--stiffness',   dest='vector', action='extend', type='string', \
                                        help='heading of column containing C11 (followed by C12, C44) field values', \
                                        metavar='<label>')
parser.add_option('-d','--direction', \
                       '--hkl',         dest='hkl', action='store', type='int', nargs=3, \
                                        help='direction of elastic modulus %default')

parser.set_defaults(vector = [])
parser.set_defaults(hkl = [1,1,1])

(options,filenames) = parser.parse_args()

if len(options.vector)== 0:
  parser.error('no data column specified...')

datainfo = {                                                               # list of requested labels per datatype
             'vector':     {'len':3,
                            'label':[]},
           }


if options.vector  != None:    datainfo['vector']['label']  += options.vector


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
  if file['name'] != 'STDIN': file['croak'].write(file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info
  table.info_append(string.replace('$Id$','\n','\\n') + \
                    '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out columns to process
  active = {}
  column = {}
  head = []

  for datatype,info in datainfo.items():
    for label in info['label']:
      foundIt = False
      for key in ['1_'+label,label]:
        if key in table.labels:
          foundIt = True
          if datatype not in active: active[datatype] = []
          if datatype not in column: column[datatype] = {}
          active[datatype].append(label)
          column[datatype][label] = table.labels.index(key)                   # remember columns of requested data
          table.labels_append('E%i%i%i(%s)'%(options.hkl[0],
                                               options.hkl[1],
                                               options.hkl[2],label)) # extend ASCII header with new labels
      if not foundIt:
        file['croak'].write('column %s not found...\n'%label)
       
# ------------------------------------------ assemble header ---------------------------------------  

  table.head_write()

# ------------------------------------------ process data ---------------------------------------  

  while table.data_read():                                                  # read next data line of ASCII table
    
    for datatype,labels in active.items():                                  # loop over vector,tensor
      for label in labels:                                                  # loop over all requested stiffnesses
        table.data_append(E_hkl(map(float,table.data[column[datatype][label]:\
                                                    column[datatype][label]+datainfo[datatype]['len']]),options.hkl))
    
    table.data_write()                                                      # output processed line

# ------------------------------------------ output result ---------------------------------------  

  table.output_flush()                                                      # just in case of buffered ASCII table

  if file['name'] != 'STDIN':
    file['input'].close()                                                   # close input ASCII table
    file['output'].close()                                                  # close output ASCII table
    os.rename(file['name']+'_tmp',file['name'])                             # overwrite old one with tmp new
