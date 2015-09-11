#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,string
import damask
from optparse import OptionParser

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [options] dfile[s]', description = """
Rename scalar, vectorial, and/or tensorial data header labels.

""", version = scriptID)

parser.add_option('-l','--label',
                  dest = 'label',
                  action = 'extend', metavar='<string LIST>',
                  help = 'column(s) to rename')
parser.add_option('-s','--substitute',
                  dest = 'substitute',
                  action = 'extend', metavar='<string LIST>',
                  help = 'new column label(s)')

parser.set_defaults(label = [],
                    substitute = [],
                   )

(options,filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except: continue
  table.report_name(scriptName,name)

# ------------------------------------------ read header ------------------------------------------  

  table.head_read()

# ------------------------------------------ process labels ---------------------------------------  

  errors  = []
  remarks = []

  if len(options.label) == 0:
    errors.append('no labels specified.')
  elif len(options.label) != len(options.substitute):
    errors.append('mismatch between number of labels ({}) and substitutes ({}).'.format(len(options.label),
                                                                                        len(options.substitute)))
  else:
    indices    = table.label_index    (options.label)
    dimensions = table.label_dimension(options.label)
    for i,index in enumerate(indices):
      if index == -1: remarks.append('label {} not present...'.format(options.label[i]))
      else:
        for j in xrange(dimensions[i]):
          table.labels[index+j] = table.labels[index+j].replace(options.label[i],options.substitute[i])

  if remarks != []: table.croak(remarks)
  if errors  != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header ---------------------------------------  

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data ------------------------------------------  

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
