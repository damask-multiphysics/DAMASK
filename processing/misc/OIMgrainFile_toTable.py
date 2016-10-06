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
Adds header to OIM grain file to make it accesible as ASCII table

""", version = scriptID)

defaults = {'labels':['1_euler','2_euler','3_euler','1_pos','2_pos','IQ','CI','Fit','GrainID']}

parser.add_option('-l', '--labels',
                  dest   = 'labels',
                  action = 'extend', metavar = '<string LIST>',
                  help   = 'lables of requested columns [{}]'.format(defaults['labels']))


(options, filenames) = parser.parse_args()

if options.labels is None:
  options.labels = defaults['labels']

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

# --- output finalization --------------------------------------------------------------------------

  table.close()                                                                                     # close ASCII table
