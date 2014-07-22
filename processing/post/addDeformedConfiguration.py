#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,string
import numpy as np
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID = '$Id$'
scriptName = scriptID.split()[1]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options file[s]', description = """
Add column(s) containing deformed configuration of requested column(s).
Operates on periodic ordered three-dimensional data sets.

""", version = string.replace(scriptID,'\n','\\n')
)

parser.add_option('-c','--coordinates', dest='coords', action='store', type='string', metavar='string', \
                                        help='column heading for coordinates [%default]')
parser.add_option('-d','--defgrad',     dest='defgrad', action='store', type='string', metavar='string', \
                                        help='heading of columns containing tensor field values')
parser.add_option('-l', '--linear',     dest='linearreconstruction', action='store_true',\
                                        help='use linear reconstruction of geometry [%default]')
parser.set_defaults(coords  = 'ip')
parser.set_defaults(defgrad = 'f' )
parser.set_defaults(linearreconstruction = False)

(options,filenames) = parser.parse_args()

datainfo = {                                                                                        # list of requested labels per datatype
             'defgrad':    {'len':9,
                            'label':[]},
           }

datainfo['defgrad']['label'].append(options.defgrad)

# ------------------------------------------ setup file handles -------------------------------------
files = []
for name in filenames:
  if os.path.exists(name):
    files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],False)                                     # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.info_append(string.replace(scriptID,'\n','\\n') + '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out dimension and resolution ----------------------------------------------
  try:
    locationCol = table.labels.index('%s.x'%options.coords)                                         # columns containing location data
  except ValueError:
    file['croak'].write('no coordinate data found...\n'%key)
    continue

  grid = [{},{},{}]
  while table.data_read():                                                                          # read next data line of ASCII table
    for j in xrange(3):
      grid[j][str(table.data[locationCol+j])] = True                                                # remember coordinate along x,y,z
  res = np.array([len(grid[0]),\
                  len(grid[1]),\
                  len(grid[2]),],'i')                                                               # resolution is number of distinct coordinates found
  geomdim = res/np.maximum(np.ones(3,'d'),res-1.0)* \
              np.array([max(map(float,grid[0].keys()))-min(map(float,grid[0].keys())),\
                        max(map(float,grid[1].keys()))-min(map(float,grid[1].keys())),\
                        max(map(float,grid[2].keys()))-min(map(float,grid[2].keys())),\
                          ],'d')                                                                    # dimension from bounding box, corrected for cell-centeredness
  if res[2] == 1:
    geomdim[2] = min(geomdim[:2]/res[:2])
  N = res.prod()

# --------------- figure out columns to process  ---------------------------------------------------
  missingColumns = False
  
  for datatype,info in datainfo.items():
    for label in info['label']:
      key = '1_%s'%label
      if key not in table.labels:
        file['croak'].write('column %s not found...\n'%key)
        missingColumns = True
      else:
        column = table.labels.index(key)                                                            # remember columns of requested data

  if missingColumns:
    continue
        
# ------------------------------------------ assemble header ---------------------------------------  
  table.labels_append(['%s_coords'%(coord+1) for coord in xrange(3)])                              # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ read deformation gradient field ----------------------- 
  table.data_rewind()
  F = np.array([0.0 for i in xrange(N*9)]).reshape([3,3]+list(res))
  idx = 0
  while table.data_read():    
    (x,y,z) = damask.util.gridLocation(idx,res)                                                     # figure out (x,y,z) position from line count
    idx += 1
    F[0:3,0:3,x,y,z] = np.array(map(float,table.data[column:column+9]),'d').reshape(3,3)

# ------------------------------------------ calculate coordinates ---------------------------------
  Favg = damask.core.math.tensorAvg(F)
  if options.linearreconstruction:
    centroids = damask.core.mesh.deformedCoordsLin(geomdim,F,Favg)
  else:
    centroids = damask.core.mesh.deformedCoordsFFT(geomdim,F,Favg)
  
# ------------------------------------------ process data ------------------------------------------
  table.data_rewind()
  idx = 0
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    (x,y,z) = damask.util.gridLocation(idx,res)                                                     # figure out (x,y,z) position from line count
    idx += 1
    table.data_append(list(centroids[:,x,y,z]))
    outputAlive = table.data_write()                                                                # output processed line
  
# ------------------------------------------ output result ---------------------------------------  
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  file['input'].close()                                                                             # close input ASCII table (works for stdin)
  file['output'].close()                                                                            # close output ASCII table (works for stdout)
  os.rename(file['name']+'_tmp',file['name'])                                                       # overwrite old one with tmp new
