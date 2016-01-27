#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options file[s]', description = """
Add column(s) containing the shape and volume mismatch resulting from given deformation gradient.
Operates on periodic three-dimensional x,y,z-ordered data sets.

""", version = scriptID)


parser.add_option('-c','--coordinates',
                  dest = 'coords',
                  type = 'string', metavar = 'string',
                  help = 'column heading of coordinates [%default]')
parser.add_option('-f','--defgrad',
                  dest = 'defgrad',
                  type = 'string', metavar = 'string ',
                  help = 'column heading of deformation gradient [%default]')
parser.add_option('--no-shape','-s',
                  dest = 'shape',
                  action = 'store_false',
                  help = 'omit shape mismatch')
parser.add_option('--no-volume','-v',
                  dest = 'volume',
                  action = 'store_false',
                  help = 'omit volume mismatch')
parser.set_defaults(coords   = 'ipinitialcoord',
                    defgrad  = 'f',
                    shape = True,
                    volume = True,
                   )

(options,filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  errors  = []
  remarks = []
  
  if table.label_dimension(options.coords) != 3:  errors.append('coordinates {} are not a vector.'.format(options.coords))
  else: colCoord = table.label_index(options.coords)

  if table.label_dimension(options.defgrad) != 9: errors.append('deformation gradient {} is not a tensor.'.format(options.defgrad))
  else: colF = table.label_index(options.defgrad)

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  if options.shape:  table.labels_append('shapeMismatch({})'.format(options.defgrad))
  if options.volume: table.labels_append('volMismatch({})'.format(options.defgrad))
  #table.head_write()

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray()

  coords = [np.unique(table.data[:,colCoord+i]) for i in xrange(3)]
  mincorner = np.array(map(min,coords))
  maxcorner = np.array(map(max,coords))
  grid   = np.array(map(len,coords),'i')
  size   = grid/np.maximum(np.ones(3,'d'), grid-1.0) * (maxcorner-mincorner)                        # size from edge to edge = dim * n/(n-1) 
  size   = np.where(grid > 1, size, min(size[grid > 1]/grid[grid > 1]))                             # spacing for grid==1 equal to smallest among other spacings

  N = grid.prod()
  
# --------------- figure out columns to process  ---------------------------------------------------
  key = '1_%s'%options.defgrad
  if key not in table.labels:
    file['croak'].write('column %s not found...\n'%key)
    continue
  else:
    column = table.labels.index(key)                                                                # remember columns of requested data

# ------------------------------------------ assemble header ---------------------------------------
  if options.shape:  table.labels_append(['shapeMismatch(%s)' %options.defgrad])
  if options.volume: table.labels_append(['volMismatch(%s)'%options.defgrad])
  table.head_write()

# ------------------------------------------ read deformation gradient field -----------------------
  table.data_rewind()
  F = np.zeros(N*9,'d').reshape([3,3]+list(grid))
  idx = 0
  while table.data_read():    
    (x,y,z) = damask.util.gridLocation(idx,grid)                                                     # figure out (x,y,z) position from line count
    idx += 1
    F[0:3,0:3,x,y,z] = np.array(map(float,table.data[column:column+9]),'d').reshape(3,3)                                               
  print 'hm'
  Favg = damask.core.math.tensorAvg(F)
  centres = damask.core.mesh.deformedCoordsFFT(size,F,Favg,[1.0,1.0,1.0])
  
  nodes   = damask.core.mesh.nodesAroundCentres(size,Favg,centres)
  if options.shape:   shapeMismatch = damask.core.mesh.shapeMismatch( size,F,nodes,centres)
  if options.volume: volumeMismatch = damask.core.mesh.volumeMismatch(size,F,nodes)

# ------------------------------------------ process data ------------------------------------------
  table.data_rewind()
  idx = 0
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    (x,y,z) = damask.util.gridLocation(idx,grid)                                                    # figure out (x,y,z) position from line count
    idx += 1
    if options.shape:  table.data_append( shapeMismatch[x,y,z])
    if options.volume: table.data_append(volumeMismatch[x,y,z])
    outputAlive = table.data_write()

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
