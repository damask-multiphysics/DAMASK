#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

def merge_dicts(*dict_args):
  """Given any number of dicts, shallow copy and merge into a new dict, with precedence going to key value pairs in latter dicts."""
  result = {}
  for dictionary in dict_args:
      result.update(dictionary)
  return result

def gradFFT(geomdim,field):
  """Calculate gradient of a vector or scalar field by transforming into Fourier space."""
  shapeFFT = np.array(np.shape(field))[0:3]
  grid     = np.array(np.shape(field)[2::-1])
  N = grid.prod()                                                                                    # field size
  n = np.array(np.shape(field)[3:]).prod()                                                           # data size

  field_fourier = np.fft.rfftn(field,axes=(0,1,2),s=shapeFFT)
  grad_fourier  = np.empty(field_fourier.shape+(3,),'c16')

  # differentiation in Fourier space
  TWOPIIMG = 2.0j*math.pi
  einsums = { 
              1:'ijkl,ijkm->ijkm',                                                                   # scalar, 1 -> 3
              3:'ijkl,ijkm->ijklm',                                                                  # vector, 3 -> 3x3
            }

  k_sk = np.where(np.arange(grid[2])>grid[2]//2,np.arange(grid[2])-grid[2],np.arange(grid[2]))/geomdim[0]
  if grid[2]%2 == 0: k_sk[grid[2]//2] = 0                                                            # Nyquist freq=0 for even grid (Johnson, MIT, 2011)

  k_sj = np.where(np.arange(grid[1])>grid[1]//2,np.arange(grid[1])-grid[1],np.arange(grid[1]))/geomdim[1]
  if grid[1]%2 == 0: k_sj[grid[1]//2] = 0                                                            # Nyquist freq=0 for even grid (Johnson, MIT, 2011)

  k_si = np.arange(grid[0]//2+1)/geomdim[2]

  kk, kj, ki = np.meshgrid(k_sk,k_sj,k_si,indexing = 'ij')
  k_s = np.concatenate((ki[:,:,:,None],kj[:,:,:,None],kk[:,:,:,None]),axis = 3).astype('c16')                           
  grad_fourier = np.einsum(einsums[n],field_fourier,k_s)*TWOPIIMG

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
parser.add_option('-l','--label',
                  dest = 'data',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'label(s) of field values')

parser.set_defaults(pos = 'pos',
                   )

(options,filenames) = parser.parse_args()

if options.data is None: parser.error('no data column specified.')

# --- define possible data types -------------------------------------------------------------------

datatypes = {
              1: {'name': 'scalar',
                  'shape': [1],
                 },
              3: {'name': 'vector',
                  'shape': [3],
                 },
            }

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
  active = []

  coordDim = table.label_dimension(options.pos)
  if coordDim != 3:
    errors.append('coordinates "{}" must be three-dimensional.'.format(options.pos))
  else: coordCol = table.label_index(options.pos)

  for me in options.data:
    dim = table.label_dimension(me)
    if dim in datatypes:
      active.append(merge_dicts({'label':me},datatypes[dim]))
      remarks.append('differentiating {} "{}"...'.format(datatypes[dim]['name'],me))
    else:
      remarks.append('skipping "{}" of dimension {}...'.format(me,dim) if dim != -1 else \
                     '"{}" not found...'.format(me) )

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  for data in active:
    table.labels_append(['{}_gradFFT({})'.format(i+1,data['label']) 
                        for i in range(coordDim*np.prod(np.array(data['shape'])))])                 # extend ASCII header with new labels
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
  for data in active:
    # we need to reverse order here, because x is fastest,ie rightmost, but leftmost in our x,y,z notation
    stack.append(gradFFT(size[::-1],
                         table.data[:,table.label_indexrange(data['label'])].
                         reshape(grid[::-1].tolist()+data['shape'])))

# ------------------------------------------ output result -----------------------------------------

  if len(stack) > 1: table.data = np.hstack(tuple(stack))
  table.data_writeArray('%.12g')

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)
