#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
from optparse import OptionParser
from collections import defaultdict
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

def gradFFT(geomdim,field):
 shapeFFT    = np.array(np.shape(field))[0:3]
 grid = np.array(np.shape(field)[2::-1])
 N = grid.prod()                                                                                    # field size
 n = np.array(np.shape(field)[3:]).prod()                                                           # data size

 if   n == 3:   dataType = 'vector'
 elif n == 1:   dataType = 'scalar'

 field_fourier = np.fft.rfftn(field,axes=(0,1,2),s=shapeFFT)
 grad_fourier  = np.empty(field_fourier.shape+(3,),'c16')

# differentiation in Fourier space
# Question: why are grid[0,1,2] normalized by geomdim[2,1,0]??
 TWOPIIMG = 2.0j*math.pi
 k_sk = np.where(np.arange(grid[2])>grid[2]//2,np.arange(grid[2])-grid[2],np.arange(grid[2]))/geomdim[0]
 if grid[2]%2 == 0: k_sk[grid[2]//2] = 0                                                            # for even grid, set Nyquist freq to 0 (Johnson, MIT, 2011)
 
 k_sj = np.where(np.arange(grid[1])>grid[1]//2,np.arange(grid[1])-grid[1],np.arange(grid[1]))/geomdim[1]
 if grid[1]%2 == 0: k_sj[grid[1]//2] = 0                                                            # for even grid, set Nyquist freq to 0 (Johnson, MIT, 2011)

 k_si = np.arange(grid[0]//2+1)/geomdim[2]
 
 kk, kj, ki = np.meshgrid(k_sk,k_sj,k_si,indexing = 'ij')
 k_s = np.concatenate((ki[:,:,:,None],kj[:,:,:,None],kk[:,:,:,None]),axis = 3).astype('c16')                           
 if dataType == 'vector':                                                                           # vector, 3 -> 3x3
   grad_fourier = np.einsum('ijkl,ijkm->ijklm',field_fourier,k_s)*TWOPIIMG
 elif dataType == 'scalar':                                                                         # scalar, 1 -> 3
   grad_fourier = np.einsum('ijkl,ijkl->ijkl',field_fourier,k_s)*TWOPIIMG

 return np.fft.irfftn(grad_fourier,axes=(0,1,2),s=shapeFFT).reshape([N,3*n])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog option(s) [ASCIItable(s)]', description = """
Add column(s) containing gradient of requested column(s).
Operates on periodic ordered three-dimensional data sets
of vector and scalar fields.

""", version = scriptID)

parser.add_option('-p','--pos','--periodiccellcenter',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'label of coordinates [%default]')
parser.add_option('-d','--data',
                  dest = 'data',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'label(s) of field values')

parser.set_defaults(pos = 'pos',
                   )

(options,filenames) = parser.parse_args()

if options.data is None: parser.error('no data column specified.')

# --- loop over input files ------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()

  remarks = []
  errors  = []
  active = defaultdict(list)

  coordDim = table.label_dimension(options.pos)
  if coordDim != 3:
    errors.append('coordinates "{}" must be three-dimensional.'.format(options.pos))
  else: coordCol = table.label_index(options.pos)

  for i,dim in enumerate(table.label_dimension(options.data)):
    me = options.data[i]
    if dim == -1:
      remarks.append('"{}" not found...'.format(me))
    elif dim == 1:
      active['scalar'].append(me)
      remarks.append('differentiating scalar "{}"...'.format(me))
    elif dim == 3:
      active['vector'].append(me)
      remarks.append('differentiating vector "{}"...'.format(me))
    else:
      remarks.append('skipping "{}" of dimension {}...'.format(me,dim))

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue


# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  for type, data in active.iteritems():
    for label in data:
      table.labels_append(['{}_gradFFT({})'.format(i+1,label) for i in range(3*table.label_dimension(label))])     # extend ASCII header with new labels
  table.head_write()

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray()

  coords = [np.unique(table.data[:,coordCol+i]) for i in range(3)]
  mincorner = np.array(map(min,coords))
  maxcorner = np.array(map(max,coords))
  grid   = np.array(map(len,coords),'i')
  size   = grid/np.maximum(np.ones(3,'d'), grid-1.0) * (maxcorner-mincorner)                        # size from edge to edge = dim * n/(n-1) 
  size   = np.where(grid > 1, size, min(size[grid > 1]/grid[grid > 1]))                             # spacing for grid==1 equal to smallest among other ones

# ------------------------------------------ process value field -----------------------------------

  stack = [table.data]
  for type, data in active.iteritems():
    for i,label in enumerate(data):
      # we need to reverse order here, because x is fastest,ie rightmost, but leftmost in our x,y,z notation
      stack.append(gradFFT(size[::-1],
                           table.data[:,table.label_indexrange(label)].
                           reshape(grid[::-1].tolist()+[table.label_dimension(label)])))

# ------------------------------------------ output result -----------------------------------------

  if len(stack) > 1: table.data = np.hstack(tuple(stack))
  table.data_writeArray('%.12g')

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)
