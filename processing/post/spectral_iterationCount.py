#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP

# -----------------------------
class extendedOption(Option):
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



# ----------------------- MAIN -------------------------------

parser = OptionParser(option_class=extendedOption, usage='%prog [options] spectralOut[s]', description = """
Generate datafile of converged iteration per increment.

""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-m','--memory',      dest='memory', action='store_true', \
                                        help='load complete file into memory [%default]')

parser.set_defaults(memory = False)

(options, filenames) = parser.parse_args()

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    if os.path.exists(name):
      (head,tail) = os.path.split(name)
      files.append({'name':name, 'input':open(name), 
                                'output':open(os.path.join(head,'iterationCount_'+os.path.splitext(tail)[0]+'.txt'), 'w')})


# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  print file['name']

# ------------------------------------------ assemble header ---------------------------------------  

  output = '1\theader\n' + \
           '\t'.join(['loadcase','increment','iterations']) + '\n'

  if options.memory:
    data = file['input'].readlines()
  else:
    data = []
    file['output'].write(output)
    output = ''
  
# ------------------------------------------ read file ---------------------------------------  

  pattern = re.compile('Loadcase\s+(\d+)\s+Increment\s+(\d+)/\d+\s+@\s+Iteration\s+(\d+)/\d+.*')
  lastIteration = 0

  for line in {True  : data,
               False : file['input']}[options.memory]:

    m = re.match(pattern, line)
    if m:
      thisLoadcase = int(m.group(1))
      thisIncrement = int(m.group(2))
      thisIteration = int(m.group(3))
      if thisIteration <= lastIteration:                                          # indicator for new increment or loadcase
        output += '\t'.join(map(str,[lastLoadcase,lastIncrement,lastIteration])) + '\n'
        if not options.memory:
          file['output'].write(output)
          output = ''
      lastLoadcase = thisLoadcase
      lastIncrement = thisIncrement
      lastIteration = thisIteration
  
  output += '\t'.join(map(str,[lastLoadcase,lastIncrement,lastIteration])) + '\n' # process last iteration (for which, of course, no further iteration can be found)
  
  file['input'].close()

# ------------------------------------------ output result ---------------------------------------  

  file['output'].write(output)

  if file['name'] != 'STDIN':
    file['output'].close
