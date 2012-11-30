#!/usr/bin/env python

import os,sys,string,damask
from optparse import OptionParser


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(usage='%prog [options] [file[s]]', description = """
Show components of given ASCIItable(s).
""" + string.replace('$Id$','\n','\\n')
)


parser.add_option('-a','--head',   dest='head',   action='store_true', help='output all heading (info + labels)')
parser.add_option('-i','--info',   dest='info',   action='store_true', help='output info lines')
parser.add_option('-l','--labels', dest='labels', action='store_true', help='output labels')
parser.add_option('-d','--data',   dest='data',   action='store_true', help='output data')
parser.add_option('-s','--sequential', dest='seq', action='store_true', help='switch to label column format')

parser.set_defaults(seq = False)
(options,filenames) = parser.parse_args()


# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':sys.stdout})

# ------------------------------------------ extract labels ---------------------------------------  

for file in files:
  table = damask.ASCIItable(file['input'],file['output'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info
  if options.head or options.info:   file['output'].write('\n'.join(table.info)+'\n')
  if options.head or options.labels: file['output'].write({True:'\n',False:'\t'}[options.seq].join(table.labels)+'\n')
  if options.data:
    while table.data_read(): table.data_write()

  table.output_flush()

  try:
    file['input'].close()
  except:
    pass
  try:
    file['output'].close()
  except:
    pass
