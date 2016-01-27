#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math,string
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# definition of element-wise p-norms for matrices

def norm(which,object):

  if which == 'Abs':                                                                                # p = 1
    return sum(map(abs, object))
  elif which == 'Frobenius':                                                                        # p = 2
    return math.sqrt(sum([x*x for x in object]))
  elif which == 'Max':                                                                              # p = inf
    return max(map(abs, object))

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing norm of requested column(s) being either vectors or tensors.

""", version = scriptID)

normChoices = ['abs','frobenius','max']
parser.add_option('-n','--norm',
                  dest = 'norm',
                  type = 'choice', choices = normChoices, metavar='string',
                  help = 'type of element-wise p-norm [frobenius] {%s}'%(','.join(map(str,normChoices))))
parser.add_option('-l','--label',
                  dest = 'label',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'heading of column(s) to calculate norm of')

parser.set_defaults(norm = 'frobenius',
                   )

(options,filenames) = parser.parse_args()

if options.label == None:
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

  errors  = []
  remarks = []
  columns = []
  dims    = []
  
  for what in options.label:
    dim = table.label_dimension(what)
    if dim < 0: remarks.append('column {} not found...'.format(what))
    else:
      dims.append(dim)
      columns.append(table.label_index(what))
      table.labels_append('norm{}({})'.format(options.norm.capitalize(),what))                    # extend ASCII header with new labels

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
    for column,dim in zip(columns,dims):
      table.data_append(norm(options.norm.capitalize(),
                             map(float,table.data[column:column+dim])))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)
