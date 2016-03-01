#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#--------------------------------------------------------------------------------------------------
def deformedCoordsFFT(F,undeformed=False):
#--------------------------------------------------------------------------------------------------
 wgt = 1.0/grid.prod()
 integrator = np.array([0.+1.j,0.+1.j,0.+1.j],'c16') * size/ 2.0 / math.pi
 step = size/grid
 
 F_fourier = np.fft.rfftn(F,axes=(0,1,2))
 coords_fourier  = np.zeros(F_fourier.shape[0:4],'c16')

 if undeformed:
   Favg=np.eye(3)
 else:
   Favg=np.real(F_fourier[0,0,0,:,:])*wgt
#--------------------------------------------------------------------------------------------------
# integration in Fourier space
 k_s = np.zeros([3],'i')
 for i in xrange(grid[2]):
   k_s[2] = i
   if(i > grid[2]//2 ): k_s[2] = k_s[2] - grid[2]
   for j in xrange(grid[1]):
     k_s[1] = j
     if(j > grid[1]//2 ): k_s[1] = k_s[1] - grid[1]
     for k in xrange(grid[0]//2+1):
       k_s[0] = k
       for m in xrange(3):
         coords_fourier[i,j,k,m] = sum(F_fourier[i,j,k,m,0:3]*k_s*integrator)
       if (any(k_s != 0)): 
         coords_fourier[i,j,k,0:3] /= -sum(k_s*k_s)

#--------------------------------------------------------------------------------------------------
# add average to scaled fluctuation and put (0,0,0) on (0,0,0)
 coords = np.fft.irfftn(coords_fourier,F.shape[0:3],axes=(0,1,2))
 
 offset_coords = np.dot(F[0,0,0,:,:],step/2.0) - scaling*coords[0,0,0,0:3]
 for z in xrange(grid[2]):
   for y in xrange(grid[1]):
     for x in xrange(grid[0]):
       coords[z,y,x,0:3] = scaling*coords[z,y,x,0:3] \
                         + offset_coords \
                         + np.dot(Favg,step*np.array([x,y,z]))

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
  centroids = deformedCoordsFFT(table.data[:,colF:colF+9].reshape(grid[2],grid[1],grid[0],3,3),
                                                                                  options.undeformed)
# ------------------------------------------ process data ------------------------------------------
  table.data_rewind()
  for z in xrange(grid[2]):
    for y in xrange(grid[1]):
      for x in xrange(grid[0]):
        table.data_read()
        table.data_append(list(centroids[z,y,x,:]))
        table.data_write()                       

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
