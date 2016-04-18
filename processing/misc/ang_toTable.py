#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [geomfile[s]]', description = """
Convert TSL/EDAX *.ang file to ASCIItable

""", version = scriptID)

(options, filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              outname = os.path.splitext(name)[0]+'.txt' if name else name,
                              buffered = False, labeled = False)
  except: continue
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name else ''))

# --- interpret header -----------------------------------------------------------------------------

  table.head_read()

# --- read comments --------------------------------------------------------------------------------

  table.info_clear()
  while table.data_read(advance = False) and table.line.startswith('#'):                            # cautiously (non-progressing) read header
    table.info_append(table.line)                                                                   # add comment to info part
    table.data_read()                                                                               # wind forward

  table.labels_clear()
  table.labels_append(['1_Euler','2_Euler','3_Euler',
                       '1_pos','2_pos',
                       'IQ','CI','PhaseID','Intensity','Fit',
                      ],                                                                            # OIM Analysis 7.2 Manual, p 403 (of 517)
                      reset = True)

# ------------------------------------------ assemble header ---------------------------------------

  table.head_write()

#--- write remainder of data file ------------------------------------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():
    outputAlive = table.data_write()

# ------------------------------------------ finalize output ---------------------------------------

  table.close()
