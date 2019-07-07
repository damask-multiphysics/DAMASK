#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser

import numpy as np

import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Append data of ASCIItable(s) row-wise.

""", version = scriptID)

parser.add_option('-a', '--add','--table',
                  dest = 'table',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'tables to add')

(options,filenames) = parser.parse_args()

if options.table is None:
  parser.error('no table specified.')


# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False)
  except: continue

  damask.util.report(scriptName,name)

  tables = []
  for addTable in options.table:
    try:    tables.append(damask.ASCIItable(name = addTable,
                                            buffered = False,
                                            readonly = True)
                         )
    except: continue

# ------------------------------------------ read headers ------------------------------------------

  table.head_read()
  for addTable in tables: addTable.head_read()

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  table.data_readArray()
  data = table.data
  for addTable in tables:
    addTable.data_readArray(table.labels(raw = True))
    data = np.vstack((data,addTable.data))
  table.data = data
  table.data_writeArray()

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
  for addTable in tables:
    addTable.close()
