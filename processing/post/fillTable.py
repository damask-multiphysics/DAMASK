#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,string,fnmatch,math,random,numpy as np
from optparse import OptionParser
import damask

scriptID = '$Id$'
scriptName = os.path.splitext(scriptID.split()[1])[0]

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
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name else ''))

# ------------------------------------------ assemble info -----------------------------------------

  table.head_read()
# ------------------------------------------ Evaluate condition ---------------------------------------  

  specials = { \
               '_row_': 0,
             }
  labels = []
  positions = []
  
  for position,label in enumerate(table.labels):
    checker = [position in table.label_indexrange(needle) for needle in options.labels]
    if any(checker) == True:
        labels.append(label)
        positions.append(position)
  interpolator = []
  condition = options.condition                                                                     # copy per file, might be altered
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
                                               table.labels.index(operand[2]))]
      except:
        parser.error('column %s not found...\n'%operand[2])
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
        table.croak('column {} not found...'.format(column))
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
    if condition == '' or eval(eval(evaluator_condition)):                                          # location of change
      for label in [x for x in options.labels if x not in set(brokenFormula)]:
        for i in xrange(len(positions)):
          table.data[positions[i]] = unravel(eval(evaluator[label]))

    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables

