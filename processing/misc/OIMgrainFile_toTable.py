#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,sys
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [file[s]]', description = """
Adds header to OIM grain file type 1 to make it accesible as ASCII table

""", version = scriptID)


(options, filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False,
                              labeled = False)
  except: continue
  damask.util.report(scriptName,name)
  table.head_read()
  data = []
  while  table.data_read():
    data.append(table.data[0:9])

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append(['1_euler','2_euler','3_euler','1_pos','2_pos','IQ','CI','Fit','GrainID'])
  table.head_write()
  for i in data:
    table.data = i
    table.data_write()

# --- output finalization --------------------------------------------------------------------------

  table.close()                                                                                     # close ASCII table
