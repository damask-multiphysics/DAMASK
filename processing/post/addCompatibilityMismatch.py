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

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options file[s]', description = """
Add column containing debug information.
Operates on periodic ordered three-dimensional data sets.

""", version = scriptID)


parser.add_option('--no-shape','-s',    dest='noShape', action='store_false',
                  help='do not calcuate shape mismatch')
parser.add_option('--no-volume','-v',   dest='noVolume', action='store_false',
                  help='do not calculate volume mismatch')
parser.add_option('-c','--coordinates', dest='coords', metavar='string',
                  help='column heading for coordinates [%default]')
parser.add_option('-f','--defgrad',     dest='defgrad', metavar='string',
                  help='column heading for coordinates [%defgrad]')
parser.set_defaults(noVolume = False)
parser.set_defaults(noShape = False)
parser.set_defaults(coords  = 'ipinitialcoord')
parser.set_defaults(defgrad = 'f')

(options,filenames) = parser.parse_args()

datainfo = {                                                                                        # list of requested labels per datatype
             'defgrad':    {'len':9,
                            'label':[]},
           }

datainfo['defgrad']['label'].append(options.defgrad)

# ------------------------------------------ setup file handles ------------------------------------
files = []
for name in filenames:
  if os.path.exists(name):
    files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

#--- loop over input files -------------------------------------------------------------------------
for file in files:
  file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],False)                                     # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out size and grid ---------------------------------------------------------
  try:
    locationCol = table.labels.index('1_%s'%options.coords)                                         # columns containing location data
  except ValueError:
    try:
      locationCol = table.labels.index('%s.x'%options.coords)                                       # columns containing location data (legacy naming scheme)
    except ValueError:
      file['croak'].write('no coordinate data (1_%s/%s.x) found...\n'%(options.coords,options.coords))
      continue

  coords = [{},{},{}]
  while table.data_read():                                                                          # read next data line of ASCII table
    for j in xrange(3):
      coords[j][str(table.data[locationCol+j])] = True                                              # remember coordinate along x,y,z
  grid = np.array([len(coords[0]),\
                   len(coords[1]),\
                   len(coords[2]),],'i')                                                            # grid is number of distinct coordinates found
  size = grid/np.maximum(np.ones(3,'d'),grid-1.0)* \
            np.array([max(map(float,coords[0].keys()))-min(map(float,coords[0].keys())),\
                      max(map(float,coords[1].keys()))-min(map(float,coords[1].keys())),\
                      max(map(float,coords[2].keys()))-min(map(float,coords[2].keys())),\
                      ],'d')                                                                        # size from bounding box, corrected for cell-centeredness

  for i, points in enumerate(grid):
    if points == 1:
      options.packing[i] = 1
      options.shift[i]   = 0
      mask = np.ones(3,dtype=bool)
      mask[i]=0
      size[i] = min(size[mask]/grid[mask])                                                          # third spacing equal to smaller of other spacing
  
  N = grid.prod()

# --------------- figure out columns to process  ---------------------------------------------------
  key = '1_%s'%datainfo['defgrad']['label'][0]
  if key not in table.labels:
    file['croak'].write('column %s not found...\n'%key)
    continue
  else:
    column = table.labels.index(key)                                                                # remember columns of requested data

# ------------------------------------------ assemble header ---------------------------------------
  if not options.noShape:  table.labels_append(['shapeMismatch(%s)' %options.defgrad])
  if not options.noVolume: table.labels_append(['volMismatch(%s)'%options.defgrad])
  table.head_write()

# ------------------------------------------ read deformation gradient field -----------------------
  table.data_rewind()
  F = np.zeros(N*9,'d').reshape([3,3]+list(grid))
  idx = 0
  while table.data_read():    
    (x,y,z) = damask.util.gridLocation(idx,grid)                                                     # figure out (x,y,z) position from line count
    idx += 1
    F[0:3,0:3,x,y,z] = np.array(map(float,table.data[column:column+9]),'d').reshape(3,3)                                               
  
  Favg = damask.core.math.tensorAvg(F)
  centres = damask.core.mesh.deformedCoordsFFT(size,F,Favg,[1.0,1.0,1.0])
  
  nodes   = damask.core.mesh.nodesAroundCentres(size,Favg,centres)
  if not options.noShape:   shapeMismatch = damask.core.mesh.shapeMismatch( size,F,nodes,centres)
  if not options.noVolume: volumeMismatch = damask.core.mesh.volumeMismatch(size,F,nodes)

# ------------------------------------------ process data ------------------------------------------
  table.data_rewind()
  idx = 0
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    (x,y,z) = damask.util.gridLocation(idx,grid)                                                    # figure out (x,y,z) position from line count
    idx += 1
    if not options.noShape:  table.data_append( shapeMismatch[x,y,z])
    if not options.noVolume: table.data_append(volumeMismatch[x,y,z])
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result ---------------------------------------  
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table
  table.output_close()                                                                              # close output ASCII table
  os.rename(file['name']+'_tmp',file['name'])                                                       # overwrite old one with tmp new
