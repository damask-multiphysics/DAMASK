#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
from optparse import OptionParser
import damask

scriptID = '$Id$'
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(usage='%prog [options] [file[s]]', description = """
Show components of given ASCIItable(s).

""", version = scriptID)


parser.add_option('-a','--head',   dest='head',   action='store_true', help='output all heading (info + labels)')
parser.add_option('-i','--info',   dest='info',   action='store_true', help='output info lines')
parser.add_option('-l','--labels', dest='labels', action='store_true', help='output labels')
parser.add_option('-d','--data',   dest='data',   action='store_true', help='output data')
parser.add_option('-c','--column', dest='col',    action='store_true', help='switch to label column format')
parser.add_option('--nolabels',    dest='nolabels',    action='store_true', help='table has no labels')

parser.set_defaults(col = False)
parser.set_defaults(nolabels = False)
(options,filenames) = parser.parse_args()


# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':sys.stdout, 'croak':sys.stderr})

# ------------------------------------------ extract labels ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  table = damask.ASCIItable(file['input'],file['output'],buffered=False,labels=not options.nolabels)        # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info
  if options.head or options.info:   file['output'].write('\n'.join(table.info)+'\n')
  if options.head or options.labels: file['output'].write({True:'\n',False:'\t'}[options.col].join(table.labels)+'\n')

# ------------------------------------------ output data ---------------------------------------  

  outputAlive = options.data
  while outputAlive and table.data_read():                                  # read next data line of ASCII table
    outputAlive = table.data_write()                                        # output line

  outputAlive and table.output_flush()

  if file['name'] != 'STDIN':
    table.input_close()  
