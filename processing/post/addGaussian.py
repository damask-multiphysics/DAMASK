#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys
import numpy as np
from optparse import OptionParser
from scipy import ndimage
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog option(s) [ASCIItable(s)]', description = """
Add column(s) containing Gaussian filtered values of requested column(s).
Operates on periodic and non-periodic ordered three-dimensional data sets.
For details see scipy.ndimage documentation.

""", version = scriptID)

parser.add_option('-p','--pos','--periodiccellcenter',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'label of coordinates [%default]')
parser.add_option('-s','--scalar',
                  dest = 'scalar',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'label(s) of scalar field values')
parser.add_option('-o','--order',
                  dest = 'order',
                  type = int,
                  metavar = 'int',
                  help = 'order of the filter')
parser.add_option('--sigma',
                  dest = 'sigma',
                  type = float,
                  metavar = 'float',
                  help = 'standard deviation')
parser.add_option('--periodic',
                  dest = 'periodic',
                  action = 'store_true',
                  help = 'assume periodic grain structure')



parser.set_defaults(pos = 'pos',
                    order = 0,
                    sigma = 1,
                    periodic = False,
                   )

(options,filenames) = parser.parse_args()

if options.scalar is None:
  parser.error('no data column specified.')

# --- loop over input files ------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  items = {
            'scalar': {'dim': 1, 'shape': [1], 'labels':options.scalar, 'active':[], 'column': []},
          }
  errors  = []
  remarks = []
  column = {}
  
  if table.label_dimension(options.pos) != 3: errors.append('coordinates {} are not a vector.'.format(options.pos))
  else: colCoord = table.label_index(options.pos)

  for type, data in items.items():
    for what in (data['labels'] if data['labels'] is not None else []):
      dim = table.label_dimension(what)
      if dim != data['dim']: remarks.append('column {} is not a {}.'.format(what,type))
      else:
        items[type]['active'].append(what)
        items[type]['column'].append(table.label_index(what))

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  for type, data in items.items():
    for label in data['active']:
      table.labels_append(['Gauss{}({})'.format(options.sigma,label)])                              # extend ASCII header with new labels
  table.head_write()

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray()

  grid,size = damask.util.coordGridAndSize(table.data[:,table.label_indexrange(options.pos)])

# ------------------------------------------ process value field -----------------------------------

  stack = [table.data]
  for type, data in items.items():
    for i,label in enumerate(data['active']):
      stack.append(ndimage.filters.gaussian_filter(table.data[:,data['column'][i]],
                                                   options.sigma,options.order,
                                                   mode = 'wrap' if options.periodic else 'nearest'
                                                   ).reshape([table.data.shape[0],1])
                   )

# ------------------------------------------ output result -----------------------------------------
  if len(stack) > 1: table.data = np.hstack(tuple(stack))
  table.data_writeArray('%.12g')

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)
