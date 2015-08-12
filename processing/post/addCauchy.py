#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing Cauchy stress based on given column(s) of deformation gradient and first Piola--Kirchhoff stress.

""", version = scriptID)

parser.add_option('-f','--defgrad',
                  dest = 'defgrad',
                  type = 'string', metavar = 'string',
                  help = 'heading of columns containing deformation gradient [%default]')
parser.add_option('-p','--stress',
                  dest = 'stress',
                  type = 'string', metavar = 'string',
                  help = 'heading of columns containing first Piola--Kirchhoff stress [%default]')

parser.set_defaults(defgrad = 'f',
                    stress  = 'p',
                   )

(options,filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name, buffered = False)
  except:
    continue
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name else ''))

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  errors = []
  column = {}
  
  for tensor in [options.defgrad,options.stress]:
    dim = table.label_dimension(tensor)
    if   dim <  0: errors.append('column {} not found.'.format(tensor))
    elif dim != 9: errors.append('column {} is not a tensor.'.format(tensor))
    else:
      column[tensor] = table.label_index(tensor)

  if errors != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append(['%i_Cauchy'%(i+1) for i in xrange(9)])                                       # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    F = np.array(map(float,table.data[column[options.defgrad]:column[options.defgrad]+9]),'d').reshape(3,3)
    P = np.array(map(float,table.data[column[options.stress ]:column[options.stress ]+9]),'d').reshape(3,3)
    table.data_append(list(1.0/np.linalg.det(F)*np.dot(P,F.T).reshape(9)))                          # [Cauchy] = (1/det(F)) * [P].[F_transpose]
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)
