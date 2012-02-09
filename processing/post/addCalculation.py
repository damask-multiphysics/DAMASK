#!/usr/bin/env python

import os,re,sys,math,string,damask
from optparse import OptionParser, Option

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

Example: distance to IP coordinates -- "math.sqrt( #ip.x#**2 + #ip.y#**2 + #ip.z#**2 )"
""" + string.replace('$Id$','\n','\\n')
)


parser.add_option('-l','--label',   dest='labels', action='extend', type='string', \
                                    help='(list of) new column labels')
parser.add_option('-f','--formula', dest='formulas', action='extend', type='string', \
                                    help='(list of) formulas corresponding to labels')

parser.set_defaults(labels= [])
parser.set_defaults(formulas= [])

(options,filenames) = parser.parse_args()

if len(options.labels) != len(options.formulas):
  parser.error('number of labels (%i) and formulas (%i) do not match'%(len(options.labels),len(options.formulas)))

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w')})

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': print file['name']

  table = damask.ASCIItable(file['input'],file['output'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info
  table.info_append(string.replace('$Id$','\n','\\n') + \
                    '\t' + ' '.join(sys.argv[1:]))

  column = {}
  evaluator = {}
  for label,formula in zip(options.labels,options.formulas):
    table.labels_append(label)
    interpolator = []
    operands = re.findall(r'#(.+?)#',formula)
    for operand in operands:
      if not operand in column:
        try:
          column[operand] = table.labels.index(operand)
        except:
          parser.error('column %s not found...\n'%operand)
      interpolator += ['float(table.data[%i])'%column[operand]]
    for operand in operands:
      formula = formula.replace('#'+operand+'#','%e')
    evaluator[label] = "'" + formula + "'%(" + ','.join(interpolator) + ")"

# ------------------------------------------ assemble header ---------------------------------------  

  table.head_write()

# ------------------------------------------ process data ---------------------------------------  

  while table.data_read():                                                  # read next data line of ASCII table
  
    for label in options.labels:
      table.data_append(eval(eval(evaluator[label])))

    table.data_write()                                                      # output processed line

# ------------------------------------------ output result ---------------------------------------  

  table.output_flush()                                                      # just in case of buffered ASCII table

  file['input'].close()                                                     # close input ASCII table
  if file['name'] != 'STDIN':
    file['output'].close                                                    # close output ASCII table
    os.rename(file['name']+'_tmp',file['name'])                             # overwrite old one with tmp new
