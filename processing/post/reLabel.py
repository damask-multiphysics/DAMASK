#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,string
import damask
from optparse import OptionParser

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [options] dfile[s]', description = """
Rename scalar, vectorial, and/or tensorial data header labels.

""", version = scriptID)

parser.add_option('-l','--label',       dest='label', action='extend', metavar='<string LIST>',
                                        help='column(s) to rename')
parser.add_option('-s','--substitute',  dest='substitute', action='extend', metavar='<string LIST>',
                                        help='new column label')

parser.set_defaults(label = [])
parser.set_defaults(substitute = [])

(options,filenames) = parser.parse_args()

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
  table.head_read()
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))                                                                                # read ASCII header info

# ------------------------------------------ process labels ---------------------------------------  

  errors = []

  if options.label == []:
    errors.append('no labels specified...')
  elif len(options.label) != len(options.substitute):
    errors.append('mismatch between number of labels ({0}) and substitutes ({1})...'.format(len(options.label),
                                                                                            len(options.substitute)))
  else:                                                                                             # tag individual candidates
    indices    = table.label_index    (options.label)
    dimensions = table.label_dimension(options.label)
    for i,index in enumerate(indices):
      if index == -1:
        errors.append('label %s not present...\n'%options.label[i])
      else:
        for j in xrange(dimensions[i]):
          table.labels[index+j] = table.labels[index+j].replace(options.label[i],options.substitute[i])

  if errors != []:
    file['croak'].write('\n'.join(errors)+'\n')
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header ---------------------------------------  

  table.head_write()

# ------------------------------------------ process data ---------------------------------------  

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result ---------------------------------------  

  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.close()                                                                                     # close ASCII tables
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
