#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

def determinant(m):
  return  +m[0]*m[4]*m[8] \
          +m[1]*m[5]*m[6] \
          +m[2]*m[3]*m[7] \
          -m[2]*m[4]*m[6] \
          -m[1]*m[3]*m[8] \
          -m[0]*m[5]*m[7]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing determinant of requested tensor column(s).

""", version = scriptID)

parser.add_option('-t','--tensor',
                  dest = 'tensor',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'heading of columns containing tensor field values')

(options,filenames) = parser.parse_args()

if options.tensor is None:
  parser.error('no data column specified.')

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  items = {
            'tensor': {'dim': 9, 'shape': [3,3], 'labels':options.tensor, 'column': []},
          }
  errors  = []
  remarks = []
  
  for type, data in items.iteritems():
    for what in data['labels']:
      dim = table.label_dimension(what)
      if dim != data['dim']: remarks.append('column {} is not a {}...'.format(what,type))
      else:
        items[type]['column'].append(table.label_index(what))
        table.labels_append('det({})'.format(what))                                                 # extend ASCII header with new labels

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for type, data in items.iteritems():
      for column in data['column']:
        table.data_append(determinant(map(float,table.data[column:
                                                           column+data['dim']])))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)