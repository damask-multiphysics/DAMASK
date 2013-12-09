#!/usr/bin/env python

import os,re,sys,math,string,damask,numpy
from optparse import OptionParser, Option

scriptID = '$Id$'
scriptName = scriptID.split()[1]

def unravel(item):
  if hasattr(item,'__contains__'): return ' '.join(map(unravel,item))
  else: return str(item)

# -----------------------------
class extendableOption(Option):
# -----------------------------
# used for definition of new option parser action 'extend', which enables to take multiple option arguments
# taken from online tutorial http://docs.python.org/library/optparse.html
  
  ACTIONS = Option.ACTIONS + ("extend",)
  STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
  TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
  ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

  def take_action(self, action, dest, opt, value, values, parser):
    if action == "extend":
      lvalue = value.split(",")
      values.ensure_value(dest, []).extend(lvalue)
    else:
      Option.take_action(self, action, dest, opt, value, values, parser)



# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) with derived values according to user defined arithmetic operation between column(s).
Columns can be specified either by label or index. Use ';' for ',' in functions.

Example: distance to IP coordinates -- "math.sqrt( #ip.x#**2 + #ip.y#**2 + round(#ip.z#;3)**2 )"
""" + string.replace(scriptID,'\n','\\n')
)


parser.add_option('-l','--label',    dest='labels', action='extend', type='string', \
                                     help='(list of) new column labels', metavar='<LIST>')
parser.add_option('-f','--formula',  dest='formulas', action='extend', type='string', \
                                     help='(list of) formulas corresponding to labels', metavar='<LIST>')
parser.set_defaults(labels= [])
parser.set_defaults(formulas= [])

(options,filenames) = parser.parse_args()

if len(options.labels) != len(options.formulas):
  parser.error('number of labels (%i) and formulas (%i) do not match'%(len(options.labels),len(options.formulas)))
  
for i in xrange(len(options.formulas)):
  options.formulas[i]=options.formulas[i].replace(';',',')

# ------------------------------------------ setup file handles ---------------------------------------  

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
  table.info_append(string.replace(scriptID,'\n','\\n') + \
                    '\t' + ' '.join(sys.argv[1:]))

  evaluator = {}
  brokenFormula = {}
  
  for label,formula in zip(options.labels,options.formulas):
    interpolator = []
    for position,column in enumerate(set(re.findall(r'#(.+?)#',formula))):                          # loop over unique set of column labels in formula
      formula = formula.replace('#'+column+'#','{%i}'%position)
      if column in specials:
        interpolator += ['specials["%s"]'%column]
      elif column.isdigit():
        if len(table.labels) > int(column):
          interpolator += ['float(table.data[%i])'%(int(column))]
        else:
          file['croak'].write('column %s not found...\n'%column)
          brokenFormula{label} = True
      else:
        try:
          interpolator += ['float(table.data[%i])'%table.labels.index(column)]
        except:
          file['croak'].write('column %s not found...\n'%column)
          brokenFormula{label} = True

    if label not in brokenFormula:  
      evaluator[label] = "'" + formula + "'.format(" + ','.join(interpolator) + ")"

# ------------------------------------------ calculate one result to get length of labels  ------
  table.data_read()
  labelLen = {}
  for label in options.labels:
    labelLen[label] = numpy.size(eval(eval(evaluator[label])))

# ------------------------------------------ assemble header ---------------------------------------  
  for label,formula in zip(options.labels,options.formulas):
    if labelLen[label] == 0:
      brokenFormula[label] = True
    if label not in brokenFormula:
      if labelLen[label] == 1:
        table.labels_append(label)
      else:
        table.labels_append(['%i_%s'%(i+1,label) for i in xrange(labelLen[label])])

  table.head_write()

# ------------------------------------------ process data ---------------------------------------  

  outputAlive = True
  table.data_rewind()

  while outputAlive and table.data_read():                                  # read next data line of ASCII table

    specials['_row_'] += 1                                                  # count row
    for label in options.labels: table.data_append(unravel(eval(eval(evaluator[label]))))
    outputAlive = table.data_write()                                        # output processed line

# ------------------------------------------ output result ---------------------------------------  

  outputAlive and table.output_flush()                                      # just in case of buffered ASCII table

  file['input'].close()                                                     # close input ASCII table
  if file['name'] != 'STDIN':
    file['output'].close()                                                  # close output ASCII table
    os.rename(file['name']+'_tmp',file['name'])                             # overwrite old one with tmp new
