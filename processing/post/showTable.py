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


parser.add_option('-d','--data',
                  dest   = 'data',
                  action = 'store_true',
                  help   = 'output data')
parser.add_option('-a','--head',
                  dest   = 'head',
                  action = 'store_true',
                  help   = 'output all heading (info + labels)')
parser.add_option('-i','--info',
                  dest   = 'info',
                  action = 'store_true',
                  help   = 'output info lines')
parser.add_option('-l','--labels',
                  dest   = 'labels',
                  action = 'store_true',
                  help   = 'output labels')
parser.add_option('-c','--column',
                  dest   = 'col',
                  action = 'store_true',
                  help   = 'print labels as one column')
parser.add_option('--nolabels',
                  dest   = 'labeled',
                  action = 'store_false',
                  help   = 'table has no labels')
parser.add_option('-t','--table',
                  dest   = 'table',
                  action = 'store_true',
                  help   = 'output heading line for proper ASCIItable format')
parser.set_defaults(head   = False,
                    info   = False,
                    labels = False,
                    data   = False,
                    col    = False,
                    labeled = True,
                    table  = False,
                   )

(options,filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = ['STDIN']

for name in filenames:
  if not (name == 'STDIN' or os.path.exists(name)): continue
  table = damask.ASCIItable(name = name, outname = None,
                            buffered = False, labeled = options.labeled, readonly = True)
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name != 'STDIN' else ''))

# ------------------------------------------ output head ---------------------------------------  

  table.head_read()
  if not (options.head or options.info):                         table.info_clear()
  if not (options.head or (options.labels and options.labeled)): table.labels_clear()

  table.head_write(header = options.table)

# ------------------------------------------ output data ---------------------------------------  

  outputAlive = options.data
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    outputAlive = table.data_write()                                                                # output line

  table.close()
