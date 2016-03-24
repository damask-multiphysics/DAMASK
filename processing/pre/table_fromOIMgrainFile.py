#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Unpack geometry files containing ranges "a to b" and/or "n of x" multiples (exclusively in one line).

""", version = scriptID)

parser.add_option('-l', '--labels',
                  dest   = 'labels',
                  help   = 'output geom file with one-dimensional data arrangement')

parser.set_defaults(labels = ['1_eulerangles','1_eulerangles','1_eulerangles',
                              '1_pos','2_pos', 'IQ', 'CI', 'Fit', 'GrainID',],
                   )

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
    data.append(table.data[0:len(options.labels)])
  
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append(options.labels)
  table.head_write()
  for i in data:
    table.data = i
    table.data_write()
