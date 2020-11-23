#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import numpy as np

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


def shapeMismatch(size,F,nodes,centres):
  """
  Routine to calculate the shape mismatch.

  shape mismatch is defined as difference between the vectors from the central point to
  the corners of reconstructed (combatible) volume element and the vectors calculated by deforming
  the initial volume element with the  current deformation gradient.
  """
  sMismatch     = np.empty(F.shape[:3])

#--------------------------------------------------------------------------------------------------
# initial positions
  delta = size/grid*.5
  coordsInitial = np.vstack((delta * np.array((-1,-1,-1)),
                             delta * np.array((+1,-1,-1)),
                             delta * np.array((+1,+1,-1)),
                             delta * np.array((-1,+1,-1)),
                             delta * np.array((-1,-1,+1)),
                             delta * np.array((+1,-1,+1)),
                             delta * np.array((+1,+1,+1)),
                             delta * np.array((-1,+1,+1))))

#--------------------------------------------------------------------------------------------------
# compare deformed original and deformed positions to actual positions
  for k in range(grid[0]):
    for j in range(grid[1]):
      for i in range(grid[2]):
       sMismatch[k,j,i] = \
         + np.linalg.norm(nodes[k,  j,  i  ,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[0,0:3]))\
         + np.linalg.norm(nodes[k+1,j,  i  ,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[1,0:3]))\
         + np.linalg.norm(nodes[k+1,j+1,i  ,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[2,0:3]))\
         + np.linalg.norm(nodes[k,  j+1,i  ,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[3,0:3]))\
         + np.linalg.norm(nodes[k,  j,  i+1,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[4,0:3]))\
         + np.linalg.norm(nodes[k+1,j,  i+1,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[5,0:3]))\
         + np.linalg.norm(nodes[k+1,j+1,i+1,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[6,0:3]))\
         + np.linalg.norm(nodes[k  ,j+1,i+1,0:3] - centres[k,j,i,0:3] - np.dot(F[k,j,i,:,:], coordsInitial[7,0:3]))
  return sMismatch


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(usage='%prog options [ASCIItable(s)]', description = """
Add column(s) containing the shape and volume mismatch resulting from given deformation gradient.
Operates on periodic three-dimensional x,y,z-ordered data sets.

""", version = scriptID)


parser.add_option('-c','--coordinates',
                  dest = 'pos',
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
parser.set_defaults(pos     = 'pos',
                    defgrad = 'f',
                    shape   = True,
                    volume  = True,
                   )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]


for name in filenames:
  damask.util.report(scriptName,name)

  table = damask.Table.load(StringIO(''.join(sys.stdin.read())) if name is None else name)
  grid,size,origin = damask.grid_filters.cell_coord0_gridSizeOrigin(table.get(options.pos))

  F = table.get(options.defgrad).reshape(tuple(grid)+(-1,),order='F').reshape(tuple(grid)+(3,3))
  nodes = damask.grid_filters.node_coord(size,F)

  if options.shape:
    centers = damask.grid_filters.cell_coord(size,F)
    shapeMismatch = shapeMismatch(size,F,nodes,centers)
    table = table.add('shapeMismatch(({}))'.format(options.defgrad),
                      shapeMismatch.reshape(-1,1,order='F'),
                      scriptID+' '+' '.join(sys.argv[1:]))

  table.save((sys.stdout if name is None else name))
