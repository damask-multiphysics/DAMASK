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
#

import os,sys,math,string,numpy as np
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id: addQuaternions.py 3828M 2015-05-15 07:28:21Z (local) $','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add Quaternions based on Crystal Frame Coordinates.

""", version = scriptID)

parser.add_option('-f','--frame',       dest='frame', nargs=4, type='string', metavar='<string string string string>',
                                        help='heading of columns containing b* vector components and three frame vectors in that order')
parser.add_option('-s','--symmetry', dest='crysym', nargs=1,type='string',metavar='<string>',
                                     help='crystal symmetry definition')                                        
parser.set_defaults(frame = None)

(options,filenames) = parser.parse_args()

if options.frame == None:
  parser.error('no data column specified...')

datainfo = {'len':4,
            'label':[]
           }

if options.frame  != None:    datainfo['label']  += options.frame

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except: continue
  table.report_name(scriptName,name)

  table.head_read()                                                                                 # read ASCII header info

# --------------- figure out columns to process  ---------------------------------------------------
  active = []
  column = {}

  for label in datainfo['label']:
    key = '1_'+label if datainfo['len'] > 1 else label                                              # non-special labels have to start with '1_'
    if key in table.labels:
        active.append(label)
        column[label] = table.labels.index(key)                                                     # remember columns of requested data
    else:
      table.croak('column %s not found...'%label)

# ------------------------------------------ assemble header ---------------------------------------
  
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append(['Q_%i'%(i+1) for i in xrange(4)])                              # extend ASCII header with new labels [1 real, 3 imaginary components]
  table.head_write()
  
# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    vec    = np.zeros([4,3])
    for i,label in enumerate(active):
      vec[i,:] = np.array(table.data[column[label]:
                                     column[label]+3])
      
    if sys.argv[1:][6]=='hexagonal':                                                      # Ensure Input matrix is orthogonal
      M=np.dot(vec[0,:],vec[2,:])
      vec[1,:]=vec[1,:]/np.linalg.norm(vec[1,:])
      vec[2,:]=M*(vec[0,:]/np.linalg.norm(vec[0,:]))
      vec[3,:]=vec[3,:]/np.linalg.norm(vec[3,:])
    else:
      vec[1,:]=vec[1,:]/np.linalg.norm(vec[1,:])
      vec[2,:]=vec[2,:]/np.linalg.norm(vec[2,:])
      vec[3,:]=vec[3,:]/np.linalg.norm(vec[3,:])
     
    
    Ori=damask.Orientation(matrix=vec[1:,:],symmetry=sys.argv[1:][6])
   
    table.data_append(np.asarray(Ori.asQuaternion()))
    
    
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.close()                                                                                     # close ASCII tables
