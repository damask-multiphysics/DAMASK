#!/usr/bin/env python2
# -*- coding: UTF-8 no BOM -*-

import os,re,sys
import math                                                                                         # noqa
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add or alter column(s) with derived values according to user-defined arithmetic operation between column(s).
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

parser.add_option('-c','--condition',
                  dest   = 'condition', metavar='string',
                  help   = 'condition to filter rows')

parser.set_defaults(condition = None,
                   )

(options,filenames) = parser.parse_args()

if options.labels is None or options.formulas is None:
  parser.error('no formulas and/or labels specified.')
if len(options.labels) != len(options.formulas):
  parser.error('number of labels ({}) and formulas ({}) do not match.'.format(len(options.labels),len(options.formulas)))

for i in xrange(len(options.formulas)):
  options.formulas[i] = options.formulas[i].replace(';',',')

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
    output = damask.ASCIItable(name = name,
                               buffered = False)
  except:
    continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header -------------------------------------------  

  table.head_read()

# -----------------------------------------------------------------------------------------------------
  specials = { \
               '_row_': 0,
             }

# ------------------------------------------ Evaluate condition ---------------------------------------
  if options.condition:  
    interpolator = []
    condition = options.condition                                                                     # copy per file, since might be altered inline
    breaker = False
  
    for position,operand in enumerate(set(re.findall(r'#(([s]#)?(.+?))#',condition))):                # find three groups
      condition = condition.replace('#'+operand[0]+'#',
                                    {  '': '{%i}'%position,
                                     's#':'"{%i}"'%position}[operand[1]])
      if operand[2] in specials:                                                                      # special label 
        interpolator += ['specials["%s"]'%operand[2]]
      else:
        try:
          interpolator += ['%s(table.data[%i])'%({  '':'float',
                                                  's#':'str'}[operand[1]],
                                                 table.label_index(operand[2]))]                     # ccould be generalized to indexrange as array lookup
        except:
          damask.util.croak('column "{}" not found.'.format(operand[2]))
          breaker = True
        
    if breaker: continue                                                                              # found mistake in condition evaluation --> next file
  
    evaluator_condition = "'" + condition + "'.format(" + ','.join(interpolator) + ")"

  else: condition = ''

# ------------------------------------------ build formulae ----------------------------------------

  evaluator = {}
  
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
        damask.util.croak('column {} not found, skipping {}...'.format(column,label))
        options.labels.remove(label)
        break

      formula = formula.replace('#'+column+'#',replacement)

    evaluator[label] = formula

    
# ------------------------------------------ process data ------------------------------------------

  firstLine   = True
  outputAlive = True

  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    specials['_row_'] += 1                                                                          # count row
    output.data_clear()
    
# ------------------------------------------ calculate one result to get length of labels  ---------

    if firstLine:
      firstLine = False
      labelDim  = {}
      for label in [x for x in options.labels]:
        labelDim[label] = np.size(eval(evaluator[label]))
        if labelDim[label] == 0: options.labels.remove(label)

# ------------------------------------------ assemble header ---------------------------------------

      output.labels_clear()
      tabLabels = table.labels()
      for label in tabLabels:
        dim = labelDim[label] if label in options.labels \
                              else table.label_dimension(label)
        output.labels_append(['{}_{}'.format(i+1,label) for i in xrange(dim)] if dim > 1 else label)

      for label in options.labels:
        if label in tabLabels: continue
        output.labels_append(['{}_{}'.format(i+1,label) for i in xrange(labelDim[label])]
                             if labelDim[label] > 1
                             else label)

      output.info = table.info
      output.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
      output.head_write()

# ------------------------------------------ process data ------------------------------------------

    for label in output.labels():
      oldIndices = table.label_indexrange(label)
      Nold = max(1,len(oldIndices))                                                                  # Nold could be zero for new columns
      Nnew = len(output.label_indexrange(label))
      output.data_append(eval(evaluator[label]) if label in options.labels and
                                                   (condition == '' or eval(eval(evaluator_condition)))
                     else np.tile([table.data[i] for i in oldIndices]
                                  if label in tabLabels
                                  else np.nan,
                                  np.ceil(float(Nnew)/Nold))[:Nnew])                                 # spread formula result into given number of columns

    outputAlive = output.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------

  table.input_close()                                                                                # close ASCII tables
  output.close()                                                                                     # close ASCII tables
  
