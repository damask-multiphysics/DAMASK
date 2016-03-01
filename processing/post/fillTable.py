#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys
import math                                                                                         # noqa
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

def unravel(item):
  if hasattr(item,'__contains__'): return ' '.join(map(unravel,item))
  else: return str(item)

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Update existing value(s) to expression(s) respecting condition.

Examples:
Replace the 2nd value of column x by value "z - y" --> fillTable -l x -c "#_row_#==2" -f "#z#-#y#" 

""", version = scriptID)

parser.add_option('-l','--label',
                  dest = 'labels',
                  action = 'extend', metavar = '<string LIST>',
                  help = '(list) of columns to be filled with formula(e)')
parser.add_option('-f','--formula',
                  dest = 'formulae',
                  action = 'extend', metavar = '<string LIST>',
                  help = '(list of) formulae corresponding to labels')
parser.add_option('-c','--condition',
                  dest   = 'condition', metavar='string',
                  help   = 'condition to filter rows')

parser.set_defaults(condition = '',
                   )

(options,filenames) = parser.parse_args()

if options.labels == None or options.formulae == None:
  parser.error('no formulae specified.')
if len(options.labels) != len(options.formulae):
  parser.error('number of labels ({}) and formulae ({}) do not match.'.format(len(options.labels),len(options.formulae)))

for i in xrange(len(options.formulae)):
  options.formulae[i] = options.formulae[i].replace(';',',')

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ assemble info -----------------------------------------

  table.head_read()
# ------------------------------------------ Evaluate condition ---------------------------------------  

  specials = { \
               '_row_': 0,
             }
  
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
                                               table.labels.index(operand[2]))]                     # ccould be generalized to indexrange as array lookup
      except:
        damask.util.croak('column %s not found.'%operand[2])
        breaker = True
        
  if breaker: continue                                                                              # found mistake in condition evaluation --> next file
  
  evaluator_condition = "'" + condition + "'.format(" + ','.join(interpolator) + ")"

#----------------------------------- Formula -------------------------------------------------------

  evaluator = {}
  brokenFormula = {}
  for label,formula in zip(options.labels,options.formulae):
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

# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))                                        # read ASCII header info
  table.labels                                                                                       # update with label set
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    specials['_row_'] += 1
    if condition == '' or eval(eval(evaluator_condition)):                                          # test row for condition
      for label in [x for x in options.labels if x not in set(brokenFormula)]:
        indices = table.label_indexrange(label)                                                     # affected columns
        probe = np.array(eval(evaluator[label]))                                                    # get formula result (scalar, array, etc.)
        container = np.tile(probe,np.ceil(float(len(indices))/probe.size))[:len(indices)]           # spread formula result into given number of columns
        for i,ind in enumerate(indices):                                                            # copy one by one as table.data is NOT a numpy array
          table.data[ind] = container[i]

    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables

