#!/usr/bin/env python

import os,re,sys,math,numpy,string,damask
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
Add column(s) containing Cauchy stress based on given column(s) of
deformation gradient and first Piola--Kirchhoff stress.

""" + string.replace(scriptID,'\n','\\n')
)


parser.add_option('-f','--defgrad',     dest='defgrad', type='string', \
                                        help='heading of columns containing deformation gradient [%default]')
parser.add_option('-p','--stress',      dest='stress', type='string', \
                                        help='heading of columns containing first Piola--Kirchhoff stress [%default]')

parser.set_defaults(defgrad = 'f')
parser.set_defaults(stress = 'p')

(options,filenames) = parser.parse_args()

if options.defgrad == None or options.stress == None:
  parser.error('missing data column...')

datainfo = {                                                               # list of requested labels per datatype
             'defgrad':    {'mandatory': True,
                            'len':9,
                            'label':[]},
             'stress':     {'mandatory': True,
                            'len':9,
                            'label':[]},
           }


datainfo['defgrad']['label'].append(options.defgrad)
datainfo['stress']['label'].append(options.stress)


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

  active = defaultdict(list)
  column = defalutdict(dict)
  missingColumns = False
  
  for datatype,info in datainfo.items():
    for label in info['label']:
      key = {True :'1_%s',
             False:'%s'   }[info['len']>1]%label
      if key not in table.labels:
        file['croak'].write('column %s not found...\n'%key)
        missingColumns |= info['mandatory']                                 # break if label is mandatory
      else:
        active[datatype].append(label)
        column[datatype][label] = table.labels.index(key)                   # remember columns of requested data

  if missingColumns:
    continue
    
  table.labels_append(['%i_Cauchy'%(i+1) 
                      for i in xrange(datainfo['stress']['len'])])          # extend ASCII header with new labels

# ------------------------------------------ assemble header ---------------------------------------

  table.head_write()

# ------------------------------------------ process data ---------------------------------------

  outputAlive = True
  while outoutAlive and table.data_read():                                  # read next data line of ASCII table
  
    F = numpy.array(map(float,table.data[column['defgrad'][active['defgrad'][0]]:
                                         column['defgrad'][active['defgrad'][0]]+datainfo['defgrad']['len']]),'d').reshape(3,3)
    P = numpy.array(map(float,table.data[column['stress'][active['stress'][0]]:
                                         column['stress'][active['stress'][0]]+datainfo['stress']['len']]),'d').reshape(3,3)

    table.data_append(list(1.0/numpy.linalg.det(F)*numpy.dot(P,F.T).reshape(9)))  # [Cauchy] = (1/det(F)) * [P].[F_transpose]
    outputAlive = table.data_write()                                        # output processed line

# ------------------------------------------ output result ---------------------------------------

  table.output_flush()                                                      # just in case of buffered ASCII table

  file['input'].close()                                                     # close input ASCII table
  if file['name'] != 'STDIN':
    file['output'].close                                                    # close output ASCII table
    os.rename(file['name']+'_tmp',file['name'])                             # overwrite old one with tmp new
