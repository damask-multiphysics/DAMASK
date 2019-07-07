#!/usr/bin/env python3

import os
from optparse import OptionParser

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Add info lines to ASCIItable header.

""", version = scriptID)

parser.add_option('-i',
                  '--info',
                  dest = 'info', action = 'extend', metavar = '<string LIST>',
                  help    = 'items to add')


(options,filenames) = parser.parse_args()

if options.info is None:
  parser.error('no info specified.')

# --- loop over input files ------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ assemble header ---------------------------------------

  table.head_read()
  table.info_append(options.info)
  table.head_write()

# ------------------------------------------ pass through data -------------------------------------

  outputAlive = True

  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close ASCII tables
