#!/usr/bin/env python

import os,re,sys,math,string,damask
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



# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Tag scalar, vectorial, and/or tensorial data header labels by specified suffix.

""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-l','--tag',         dest='tag', type='string', \
                                        help='tag to use as suffix for labels')
parser.add_option('-v','--vector',      dest='vector', action='extend', type='string', \
                                        help='heading of columns containing 3x1 vector field values')
parser.add_option('-t','--tensor',      dest='tensor', action='extend', type='string', \
                                        help='heading of columns containing 3x3 tensor field values')
parser.add_option('-s','--special',     dest='special', action='extend', type='string', \
                                        help='heading of columns containing field values of special dimension')
parser.add_option('-d','--dimension',   dest='N', action='store', type='int', \
                                        help='dimension of special field values [%default]')

parser.set_defaults(tag = '')
parser.set_defaults(vector = [])
parser.set_defaults(tensor = [])
parser.set_defaults(special = [])
parser.set_defaults(N = 1)

(options,filenames) = parser.parse_args()

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


# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w')})

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': print file['name']

  table = damask.ASCIItable(file['input'],file['output'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info
  table.info_append(string.replace('$Id$','\n','\\n') + \
                    '\t' + ' '.join(sys.argv[1:]))

# ------------------------------------------ process labels ---------------------------------------  

  if options.vector == [] and options.tensor == [] and options.special == []: # default to tagging all labels
    for i,label in enumerate(table.labels):
      table.labels[i] += options.tag
  else:                                                                       # tag individual candidates
    for datatype,info in datainfo.items():
      for label in info['label']:
        key = {True :'1_%s',
               False:'%s'   }[info['len']>1]%label
        if key not in table.labels:
          sys.stderr.write('column %s not found...\n'%key)
        else:
          offset = table.labels.index(key)
          for i in xrange(info['len']):
            table.labels[offset+i] += options.tag

# ------------------------------------------ assemble header ---------------------------------------  

  table.head_write()

# ------------------------------------------ process data ---------------------------------------  

  outputAlive = True
  while outputAlive and table.data_read():                                  # read next data line of ASCII table
    outputAlive = table.data_write()                                        # output processed line

# ------------------------------------------ output result ---------------------------------------  

  outputAlive and table.output_flush()                                      # just in case of buffered ASCII table

  file['input'].close()                                                     # close input ASCII table
  if file['name'] != 'STDIN':
    file['output'].close                                                    # close output ASCII table
    os.rename(file['name']+'_tmp',file['name'])                             # overwrite old one with tmp new
