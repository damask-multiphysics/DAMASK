#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add data in column(s) of mapped ASCIItable selected from the row indexed by the value in a mapping column.
Row numbers start at 1.

""", version = scriptID)

parser.add_option('--index',
                  dest = 'index',
                  type = 'string', metavar = 'string',
                  help = 'column label containing row index')
parser.add_option('-o','--offset',
                  dest = 'offset',
                  type = 'int', metavar = 'int',
                  help = 'constant offset for index column value [%default]')
parser.add_option('-l','--label',
                  dest = 'label',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'column label(s) to be appended')
parser.add_option('-a','--asciitable',
                  dest = 'asciitable',
                  type = 'string', metavar = 'string',
                  help = 'indexed ASCIItable')

parser.set_defaults(offset = 0,
                   )

(options,filenames) = parser.parse_args()

if options.label is None:
  parser.error('no data columns specified.')
if options.index is None:
  parser.error('no index column given.')

# ------------------------------------------ process indexed ASCIItable ---------------------------

if options.asciitable is not None and os.path.isfile(options.asciitable):

  indexedTable = damask.ASCIItable(name = options.asciitable,
                                   buffered = False,
                                   readonly = True) 
  indexedTable.head_read()                                                                          # read ASCII header info of indexed table
  missing_labels = indexedTable.data_readArray(options.label)
  indexedTable.close()                                                                              # close input ASCII table

  if len(missing_labels) > 0:
    damask.util.croak('column{} {} not found...'.format('s' if len(missing_labels) > 1 else '',', '.join(missing_labels)))

else:
  parser.error('no indexed ASCIItable given.')

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  errors = []

  indexColumn = table.label_index(options.index)  
  if indexColumn <  0: errors.append('index column {} not found.'.format(options.index))

  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append(indexedTable.labels(raw = True))                                              # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    try:
      table.data_append(indexedTable.data[int(round(float(table.data[indexColumn])))+options.offset-1]) # add all mapped data types
    except IndexError:
      table.data_append(np.nan*np.ones_like(indexedTable.data[0]))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
