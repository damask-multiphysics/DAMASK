#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,collections
import math                                                                                         # noqa
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

def listify(x):
  return x if isinstance(x, collections.Iterable) else [x]


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add or alter column(s) with derived values according to user-defined arithmetic operation between column(s).
Column labels are tagged by '#label#' in formulas. Use ';' for ',' in functions.
Numpy is available as 'np'.

Special variables: #_row_# -- row index
Examples:
(1) magnitude of vector -- "np.linalg.norm(#vec#)"
(2) rounded root of row number -- "round(math.sqrt(#_row_#);3)"

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

for i in range(len(options.formulas)):
  options.formulas[i] = options.formulas[i].replace(';',',')

# ------------------------------------- loop over input files --------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header -------------------------------------------  

  table.head_read()

# --------------------------------------------------------------------------------------------------
  specials = { \
               '_row_': 0,
             }

# --------------------------------------- evaluate condition ---------------------------------------
  if options.condition is not None:
    interpolator = []
    condition = options.condition                                                                   # copy per file, since might be altered inline
    breaker = False
  
    for position,operand in enumerate(set(re.findall(r'#(([s]#)?(.+?))#',condition))):              # find three groups
      condition = condition.replace('#'+operand[0]+'#',
                                    {  '': '{%i}'%position,
                                     's#':'"{%i}"'%position}[operand[1]])
      if operand[2] in specials:                                                                    # special label 
        interpolator += ['specials["%s"]'%operand[2]]
      else:
        try:
          interpolator += ['%s(table.data[%i])'%({  '':'float',
                                                  's#':'str'}[operand[1]],
                                                 table.label_index(operand[2]))]                    # could be generalized to indexrange as array lookup
        except:
          damask.util.croak('column "{}" not found.'.format(operand[2]))
          breaker = True
        
    if breaker: continue                                                                            # found mistake in condition evaluation --> next file
  
    evaluator_condition = "'" + condition + "'.format(" + ','.join(interpolator) + ")"

# ------------------------------------------ build formulas ----------------------------------------

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

# ---------------------------- separate requested labels into old and new --------------------------

  veterans = list(set(options.labels)&set(table.labels(raw=False)+table.labels(raw=True)) )         # intersection of requested and existing
  newbies  = list(set(options.labels)-set(table.labels(raw=False)+table.labels(raw=True)) )         # requested but not existing
    
# ------------------------------------------ process data ------------------------------------------

  firstLine   = True
  outputAlive = True

  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    specials['_row_'] += 1                                                                          # count row
    
    if firstLine:
      firstLine = False

# ---------------------------- line 1: determine dimension of formulas -----------------------------

      resultDim  = {}
      for label in list(options.labels):                                                            # iterate over stable copy
        resultDim[label] = np.size(eval(evaluator[label]))                                          # get dimension of formula[label]
        if resultDim[label] == 0: options.labels.remove(label)                                      # remove label if invalid result
      
      for veteran in list(veterans):
        if resultDim[veteran] != table.label_dimension(veteran):
          damask.util.croak('skipping {} due to inconsistent dimension...'.format(veteran))
          veterans.remove(veteran)                                                                  # discard culprit

# ----------------------------------- line 1: assemble header --------------------------------------

      for newby in newbies:
        table.labels_append(['{}_{}'.format(i+1,newby) for i in range(resultDim[newby])] 
                             if resultDim[newby] > 1 else newby)

      table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
      table.head_write()

# -------------------------------------- evaluate formulas -----------------------------------------

    if options.condition is None or eval(eval(evaluator_condition)):                                # condition for veteran replacement fulfilled
      for veteran in veterans:                                                                      # evaluate formulas that overwrite
        table.data[table.label_index(veteran):
                   table.label_index(veteran)+table.label_dimension(veteran)] = \
                   listify(eval(evaluator[veteran]))
    
    for newby in newbies:                                                                           # evaluate formulas that append
      table.data_append(listify(eval(evaluator[newby])))

    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------- output finalization ----------------------------------------

  table.close()                                                                                     # close ASCII table

