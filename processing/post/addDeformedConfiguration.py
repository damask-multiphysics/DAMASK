#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

#--------------------------------------------------------------------------------------------------
def deformedCoordsFFT(F,undeformed=False):
#--------------------------------------------------------------------------------------------------
 wgt = 1.0/grid.prod()
 integrator = grid[::-1] / 2.0 / math.pi
 step = size[::-1]/grid[::-1]

 F_fourier = np.fft.fftpack.rfftn(F,axes=(0,1,2))
 if undeformed:
   Favg=np.eye(3)
 else:
   Favg=np.real(F_fourier[0,0,0,:,:])*wgt
 coords_fourier  = np.zeros(F_fourier.shape[:-1],'c16')

#--------------------------------------------------------------------------------------------------
# integration in Fourier space
 k_s = np.zeros([3],'i')

 for k in xrange(grid[2]):
   k_s[2] = k
   if(k > grid[0]/2 ): k_s[2] = k_s[2] - grid[2]
   for j in xrange(grid[1]):
     k_s[1] = j
     if(j > grid[1]/2 ): k_s[1] = k_s[1] - grid[1]
     for i in xrange(grid[0]/2+1):
       k_s[0] = i
       for m in xrange(3):
         coords_fourier[k,j,i,m] = sum(F_fourier[k,j,i,0:3,m]*k_s*\
                                                 np.array([0.+1.j,0.+1.j,0.+1.j],'c16')*integrator)
       if (any(k_s != 0)): 
         coords_fourier[k,j,i,0:3] /= -sum(k_s*k_s)
#--------------------------------------------------------------------------------------------------
# add average to scaled fluctuation and put (0,0,0) on (0,0,0)
 coords = np.fft.fftpack.irfftn(coords_fourier,axes=(0,1,2))
 offset_coords = np.dot(F[0,0,0,0:3,0:3],step/2.0) - scaling*coords[0,0,0,0:3]
 for k in xrange(grid[2]):
   for j in xrange(grid[1]):
     for i in xrange(grid[0]):
       coords[k,j,i,0:3] = scaling*coords[k,j,i,0:3] \
                         + offset_coords \
                         + np.dot(Favg,step*np.array([i,j,k]))
 return coords

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options file[s]', description = """
Add deformed configuration of given initial coordinates.
Operates on periodic three-dimensional x,y,z-ordered data sets.

""", version = scriptID)

parser.add_option('-f', '--defgrad',dest='defgrad', metavar = 'string',
                                    help='heading of deformation gradient columns [%default]')
parser.add_option('--reference',    dest='undeformed', action='store_true',
                                    help='map results to reference (undeformed) average configuration [%default]')
parser.add_option('--scaling',      dest='scaling', action='extend', metavar = '<float LIST>',
                                    help='scaling of fluctuation')
parser.add_option('-u', '--unitlength', dest='unitlength', type='float', metavar = 'float',
                                    help='set unit length for 2D model [%default]')
parser.add_option('--coordinates',  dest='coords', metavar='string',
                                    help='column heading for coordinates [%default]')

parser.set_defaults(defgrad = 'f')
parser.set_defaults(coords = 'ipinitialcoord')
parser.set_defaults(scaling = [])
parser.set_defaults(undeformed = False)
parser.set_defaults(unitlength = 0.0)

(options,filenames) = parser.parse_args()

options.scaling += [1.0 for i in xrange(max(0,3-len(options.scaling)))]
scaling = map(float, options.scaling)


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

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray()

  coords = [np.unique(table.data[:,colCoord+i]) for i in xrange(3)]
  mincorner = np.array(map(min,coords))
  maxcorner = np.array(map(max,coords))
  grid   = np.array(map(len,coords),'i')
  size   = grid/np.maximum(np.ones(3,'d'), grid-1.0) * (maxcorner-mincorner)                        # size from edge to edge = dim * n/(n-1) 
  size   = np.where(grid > 1, size, min(size[grid > 1]/grid[grid > 1]))                             # spacing for grid==1 equal to smallest among other spacings

  N = grid.prod()

  if N != len(table.data): errors.append('data count {} does not match grid {}x{}x{}.'.format(N,*grid))
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue
  
# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  for coord in xrange(3):
    label = '{}_{}_{}'.format(coord+1,options.defgrad,options.coords)
    if np.any(scaling) != 1.0: label+='_{}_{}_{}'.format(scaling)
    if options.undeformed: label+='_undeformed'
    table.labels_append([label])                                                                    # extend ASCII header with new labels
  table.head_write()

# ------------------------------------------ read deformation gradient field -----------------------
  centroids = deformedCoordsFFT(table.data[:,colF:colF+9].reshape(grid[0],grid[1],grid[2],3,3),
                                                                                  options.undeformed)

# ------------------------------------------ process data ------------------------------------------
  table.data_rewind()
  idx = 0
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    (x,y,z) = damask.util.gridLocation(idx,grid)                                                    # figure out (x,y,z) position from line count
    idx += 1
    table.data_append(list(centroids[z,y,x,:]))
    outputAlive = table.data_write()                       

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
