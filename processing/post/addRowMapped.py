#!/usr/bin/env python2.7
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
Add data in column(s) of second ASCIItable selected from the row indexed by the value in a mapping column.

""", version = scriptID)

parser.add_option('-c','--map',
                  dest = 'map',
                  type = 'string', metavar = 'string',
                  help = 'column label containing row mapping')
parser.add_option('-o','--offset',
                  dest = 'offset',
                  type = 'int', metavar = 'int',
                  help = 'offset between mapping column value and actual row in mapped table [%default]')
parser.add_option('-l','--label',
                  dest = 'label',
                  action = 'extend', metavar = '<string LIST>',
                  help='column label(s) to be appended')
parser.add_option('-a','--asciitable',
                  dest = 'asciitable',
                  type = 'string', metavar = 'string',
                  help = 'indexed ASCIItable')

parser.set_defaults(offset = 0,
                   )

(options,filenames) = parser.parse_args()

if options.label is None:
  parser.error('no data columns specified.')
if options.map is None:
  parser.error('no mapping column given.')

# ------------------------------------------ process mapping ASCIItable ---------------------------

if options.asciitable is not None and os.path.isfile(options.asciitable):

  mappedTable = damask.ASCIItable(name = options.asciitable,
                                  buffered = False,
                                  readonly = True) 
  mappedTable.head_read()                                                                           # read ASCII header info of mapped table
  missing_labels = mappedTable.data_readArray(options.label)
  mappedTable.close()                                                                               # close mapped input ASCII table

  if len(missing_labels) > 0:
    damask.util.croak('column{} {} not found...'.format('s' if len(missing_labels) > 1 else '',', '.join(missing_labels)))

else:
  parser.error('no mapped ASCIItable given.')

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

  mappedColumn = table.label_index(options.map)  
  if mappedColumn <  0: errors.append('mapping column {} not found.'.format(options.map))

  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append(mappedTable.labels(raw = True))                                              # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    try:
      table.data_append(mappedTable.data[int(round(float(table.data[mappedColumn])))+options.offset-1]) # add all mapped data types
    except IndexError:
      table.data_append(np.nan*np.ones_like(mappedTable.data[0]))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
