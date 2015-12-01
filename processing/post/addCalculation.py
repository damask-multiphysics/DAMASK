#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,string
import math                                                                                         # flake8: noqa
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

def unravel(item):
  if hasattr(item,'__contains__'): return ' '.join(map(unravel,item))
  else: return str(item)

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) with derived values according to user-defined arithmetic operation between column(s).
Column labels are tagged by '#label#' in formulas. Use ';' for ',' in functions.
Numpy is available as np.

Special variables: #_row_# -- row index
Examples: (1) magnitude of vector -- "np.linalg.norm(#vec#)" (2) rounded root of row number -- "round(math.sqrt(#_row_#);3)"

""", version = scriptID)

parser.add_option('-l','--label',
                  dest = 'labels',
                  action = 'extend', metavar = '<string LIST>',
                  help = '(list of) new column labels')
parser.add_option('-f','--formula',
                  dest = 'formulas',
                  action = 'extend', metavar = '<string LIST>',
                  help = '(list of) formulas corresponding to labels')

(options,filenames) = parser.parse_args()

if options.labels == None or options.formulas == None:
  parser.error('no formulas and/or labels specified.')
if len(options.labels) != len(options.formulas):
  parser.error('number of labels ({}) and formulas ({}) do not match.'.format(len(options.labels),len(options.formulas)))

for i in xrange(len(options.formulas)):
  options.formulas[i] = options.formulas[i].replace(';',',')

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name, buffered = False)
  except:
    continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header -------------------------------------------  

  table.head_read()

# ------------------------------------------ build formulae ----------------------------------------

  specials = { \
               '_row_': 0,
             }

  evaluator = {}
  brokenFormula = {}
  
  for label,formula in zip(options.labels,options.formulas):
    for column in re.findall(r'#(.+?)#',formula):                                                   # loop over column labels in formula
      idx = table.label_index(column)
      dim = table.label_dimension(column)
      if column in specials:
        replacement = 'specials["{}"]'.format(column)
      elif dim == 1:                                                                                # scalar input
        replacement = 'float(table.data[{}])'.format(idx)                                           # take float value of data column
      elif dim > 1:                                                                                 # multidimensional input (vector, tensor, etc.)
        replacement = 'np.array(table.data[{}:{}],dtype=float)'.format(idx,idx+dim)                 # use (flat) array representation
      else:
        damask.util.croak('column {} not found...'.format(column))
        brokenFormula[label] = True
        break

      formula = formula.replace('#'+column+'#',replacement)
    
    if label not in brokenFormula:
      evaluator[label] = formula

# ------------------------------------------ process data ------------------------------------------

  firstLine   = True
  outputAlive = True

  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    specials['_row_'] += 1                                                                          # count row

# ------------------------------------------ calculate one result to get length of labels  ---------

    if firstLine:
      firstLine = False
      labelDim  = {}
      for label in [x for x in options.labels if x not in set(brokenFormula)]:
        labelDim[label] = np.size(eval(evaluator[label]))
        if labelDim[label] == 0: brokenFormula[label] = True

# ------------------------------------------ assemble header ---------------------------------------

        if label not in brokenFormula:
          table.labels_append(['{}_{}'.format(i+1,label) for i in xrange(labelDim[label])] if labelDim[label] > 1
                               else label)

      table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
      table.head_write()

# ------------------------------------------ process data ------------------------------------------

    for label in [x for x in options.labels if x not in set(brokenFormula)]:
      table.data_append(unravel(eval(evaluator[label])))

    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
