#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,string
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
Add column(s) with derived values according to user defined arithmetic operation between column(s).
Columns can be specified either by label or index. Use ';' for ',' in functions.
Numpy is available as np.

Example: distance to IP coordinates -- "math.sqrt( #ip.x#**2 + #ip.y#**2 + round(#ip.z#;3)**2 )"

""", version = scriptID)

parser.add_option('-l','--label',   dest='labels', action='extend', metavar='<string LIST>',
                                    help='(list of) new column labels')
parser.add_option('-f','--formula', dest='formulas', action='extend', metavar='<string LIST>',
                                    help='(list of) formulas corresponding to labels')
parser.set_defaults(labels= [])
parser.set_defaults(formulas= [])

(options,filenames) = parser.parse_args()

if len(options.labels) != len(options.formulas):
  parser.error('number of labels (%i) and formulas (%i) do not match'%(len(options.labels),len(options.formulas)))
  
for i in xrange(len(options.formulas)):
  options.formulas[i]=options.formulas[i].replace(';',',')

# ------------------------------------------ setup file handles ------------------------------------
files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  specials = { \
               '_row_': 0,
             }

  table = damask.ASCIItable(file['input'],file['output'],False)                                     # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

  evaluator = {}
  brokenFormula = {}
  
  for label,formula in zip(options.labels,options.formulas):
    interpolator = []
    for column in re.findall(r'#(.+?)#',formula):                                                   # loop over column labels in formula
      formula = formula.replace('#'+column+'#','%f')
      if column in specials:
        interpolator += ['specials["%s"]'%column]
      elif column.isdigit():
        if len(table.labels) > int(column):
          interpolator += ['float(table.data[%i])'%(int(column))]
        else:
          file['croak'].write('column %s not found...\n'%column)
          brokenFormula[label] = True
      else:
        try:
          interpolator += ['float(table.data[%i])'%table.labels.index(column)]
        except:
          file['croak'].write('column %s not found...\n'%column)
          brokenFormula[label] = True
    
    if label not in brokenFormula:
      evaluator[label] = "'" + formula + "'%(" + ','.join(interpolator) + ")"

# ------------------------------------------ process data ------------------------------------------
  firstLine=True 
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    specials['_row_'] += 1                                                                          # count row
# ------------------------------------------ calculate one result to get length of labels  ---------
    if firstLine:
      labelLen = {}
      for label in options.labels:
        labelLen[label] = np.size(eval(eval(evaluator[label])))

# ------------------------------------------ assemble header ---------------------------------------
      for label,formula in zip(options.labels,options.formulas):
        if labelLen[label] == 0:
          brokenFormula[label] = True
        if label not in brokenFormula:
          table.labels_append({True:['%i_%s'%(i+1,label) for i in xrange(labelLen[label])],
                               False:label}[labelLen[label]>1] )
      table.head_write()
      firstLine = False

    for label in options.labels: table.data_append(unravel(eval(eval(evaluator[label]))))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table (works for stdin)
  table.output_close()                                                                              # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
