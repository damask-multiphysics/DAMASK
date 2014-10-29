#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,string,fnmatch
from optparse import OptionParser
import damask

scriptID = '$Id$'
scriptName = scriptID.split()[1][:-3]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Filter rows according to condition and columns by either white or black listing.

Examples:
Every odd row if x coordinate is positive -- " #ip.x# >= 0.0 and #_row_#%2 == 1 ).
All rows where label 'foo' equals 'bar' -- " #foo# == \"bar\" "

""", version = scriptID)

parser.add_option('-w','--white',   dest='whitelist', action='extend', metavar='<string LIST>',
                                    help='white list of column labels (a,b,c,...)')
parser.add_option('-b','--black',   dest='blacklist', action='extend', metavar='<string LIST>',
                                    help='black list of column labels (a,b,c,...)')
parser.add_option('-c','--condition', dest='condition', metavar='string',
                                    help='condition to filter rows')

parser.set_defaults(whitelist = [])
parser.set_defaults(blacklist = [])
parser.set_defaults(condition = '')

(options,filenames) = parser.parse_args()

if filenames == []:
  filenames = ['STDIN']

#--- loop over input files -------------------------------------------------------------------------
for name in filenames:
  if name == 'STDIN':
    file = {'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m\n')
  else:
    if not os.path.exists(name): continue
    file = {'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],False)                                     # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

  specials = { \
               '_row_': 0,
             }
  labels = []
  positions = []

  for position,label in enumerate(table.labels):
    if    (options.whitelist == [] or     any([fnmatch.fnmatch(label,needle) for needle in options.whitelist])) \
      and (options.blacklist == [] or not any([fnmatch.fnmatch(label,needle) for needle in options.blacklist])):  # a label to keep?
      labels.append(label)                                                                          # remember name...
      positions.append(position)                                                                    # ...and position

  interpolator = []
  condition = options.condition                                                                     # copy per file, might be altered
  for position,operand in enumerate(set(re.findall(r'#(([s]#)?(.+?))#',condition))):                # find three groups
    condition = condition.replace('#'+operand[0]+'#',
                                          {  '': '{%i}'%position,
                                           's#':'"{%i}"'%position}[operand[1]])
    if operand[2] in specials:                                                                      # special label ?
      interpolator += ['specials["%s"]'%operand[2]]
    else:
      try:
        interpolator += ['%s(table.data[%i])'%({  '':'float',
                                                's#':'str'}[operand[1]],
                                               table.labels.index(operand[2]))]
      except:
        parser.error('column %s not found...\n'%operand[2])

  evaluator = "'" + condition + "'.format(" + ','.join(interpolator) + ")"
  
# ------------------------------------------ assemble header ---------------------------------------
  table.labels = labels                                                                             # update with new label set
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    specials['_row_'] += 1                                                                          # count row
    if condition == '' or eval(eval(evaluator)):                                                    # valid row ?
      table.data = [table.data[position] for position in positions]                                 # retain filtered columns
      outputAlive = table.data_write()                                                              # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table (works for stdin)
  table.output_close()                                                                              # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
