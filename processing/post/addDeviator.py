#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

oneThird = 1.0/3.0

def deviator(m,spherical = False):                                                                                    # Carefull, do not change the value of m (its intent(inout)!)
  sph = oneThird*(m[0]+m[4]+m[8])
  dev = [
           m[0]-sph,  m[1],     m[2],
           m[3],      m[4]-sph, m[5],
           m[6],      m[7],     m[8]-sph,
        ]
  return dev,sph if spherical else dev

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing deviator of requested tensor column(s).

""", version = scriptID)

parser.add_option('-t','--tensor',
                  dest = 'tensor',
                  action = 'extend', metavar='<string LIST>',
                  help = 'heading of columns containing tensor field values')
parser.add_option('-s','--spherical',
                  dest = 'spherical',
                  action = 'store_true',
                  help = 'report spherical part of tensor (hydrostatic component, pressure)')

(options,filenames) = parser.parse_args()

if options.tensor == None:
  parser.error('no data column specified...')

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name, buffered = False)
  except:
    continue
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name else ''))

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  items = {
            'tensor': {'dim': 9, 'shape': [3,3], 'labels':options.tensor, 'active':[], 'column': []},
          }
  errors  = []
  remarks = []
  column = {}
  
  for type, data in items.iteritems():
    for what in data['labels']:
      dim = table.label_dimension(what)
      if dim != data['dim']: remarks.append('column {} is not a {}.'.format(what,type))
      else:
        items[type]['active'].append(what)
        items[type]['column'].append(table.label_index(what))

  if remarks != []: table.croak(remarks)
  if errors  != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  for type, data in items.iteritems():
    for label in data['active']:
      table.labels_append(['{}_dev({})'.format(i+1,label) for i in xrange(data['dim'])] + \
                         (['sph({})'.format(label)] if options.spherical else []))                  # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for type, data in items.iteritems():
      for column in data['column']:
        table.data_append(deviator(map(float,table.data[column:
                                                        column+data['dim']]),options.spherical))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)
