#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-
# Copyright 2011-14 Max-Planck-Institut für Eisenforschung GmbH
#
# This file is part of DAMASK,
# the Düsseldorf Advanced Material Simulation Kit.
#
# DAMASK is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DAMASK is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with DAMASK. If not, see <http://www.gnu.org/licenses/>.


import os,sys,math,string,numpy as np
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id: addOIMTransCoord.py 3828M 2015-05-06 07:28:21Z (local) $','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Transform X,Y,Z,F APS BeamLine 34 coordinates to x,y,z APS strain coordinates.

""", version = scriptID)

parser.add_option('-f','--frame',       dest='frame', nargs=4, type='string', metavar='<string string string string>',
                                        help='APS X,Y,Z coords, and depth F')
parser.set_defaults(frame = None)

(options,filenames) = parser.parse_args()

if options.frame == None:
  parser.error('no data column specified...')


datainfo = {'len':3,
            'label':[]
           }

if options.frame  != None:    datainfo['label']  += options.frame

# --- loop over input files -------------------------------------------------------------------------
if filenames == []:
  filenames = ['STDIN']

for name in filenames:
  if name == 'STDIN':
    file = {'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m\n')
  else:
    if not os.path.exists(name): continue
    file = {'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],buffered=False)                            # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out columns to process  ---------------------------------------------------
  active = []
  column = {}
  columnMissing = False

  for label in datainfo['label']:    
    key = label
    if key in table.labels:
        active.append(label)
        column[label] = table.labels.index(key)                                                     # remember columns of requested data
    else:
      file['croak'].write('column %s not found...\n'%label)
      columnMissing = True
      
  if columnMissing: continue

# ------------------------------------------ assemble header ---------------------------------------
  table.labels_append(['%i_coord'%(i+1) for i in xrange(3)])                                       # extend ASCII header with new labels
  table.head_write()
  
# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  vec         = np.zeros(4)
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for i,label in enumerate(active):
      vec[i] = table.data[column[label]]                                  
    table.data_append(np.array([vec[0],-math.sqrt(0.5)*(vec[1]+vec[2]),-vec[3]]))               # change X,Y,Z (equivalent to x,y,z crystal orientation basis) to strain coordinate system (x radially outward of beam, y down on sample surface, z normal up from sample surface)
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table (works for stdin)
  table.output_close()                                                                              # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new

