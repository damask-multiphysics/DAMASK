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

def curlFFT(geomdim,field):
  """Calculate curl of a vector or tensor field by transforming into Fourier space."""
  shapeFFT    = np.array(np.shape(field))[0:3]
  grid = np.array(np.shape(field)[2::-1])
  N = grid.prod()                                                                                    # field size
  n = np.array(np.shape(field)[3:]).prod()                                                           # data size

  field_fourier = np.fft.rfftn(field,axes=(0,1,2),s=shapeFFT)
  curl_fourier  = np.empty(field_fourier.shape,'c16')

  # differentiation in Fourier space
  TWOPIIMG = 2.0j*math.pi
  einsums = { 
              3:'slm,ijkl,ijkm->ijks',                                                               # vector, 3 -> 3
              9:'slm,ijkl,ijknm->ijksn',                                                             # tensor, 3x3 -> 3x3
            }
  k_sk = np.where(np.arange(grid[2])>grid[2]//2,np.arange(grid[2])-grid[2],np.arange(grid[2]))/geomdim[0]
  if grid[2]%2 == 0: k_sk[grid[2]//2] = 0                                                            # Nyquist freq=0 for even grid (Johnson, MIT, 2011)

  k_sj = np.where(np.arange(grid[1])>grid[1]//2,np.arange(grid[1])-grid[1],np.arange(grid[1]))/geomdim[1]
  if grid[1]%2 == 0: k_sj[grid[1]//2] = 0                                                            # Nyquist freq=0 for even grid (Johnson, MIT, 2011)

  k_si = np.arange(grid[0]//2+1)/geomdim[2]

  kk, kj, ki = np.meshgrid(k_sk,k_sj,k_si,indexing = 'ij')
  k_s = np.concatenate((ki[:,:,:,None],kj[:,:,:,None],kk[:,:,:,None]),axis = 3).astype('c16')

  e = np.zeros((3, 3, 3))
  e[0, 1, 2] = e[1, 2, 0] = e[2, 0, 1] = 1.0                                                         # Levi-Civita symbols 
  e[0, 2, 1] = e[2, 1, 0] = e[1, 0, 2] = -1.0

  curl_fourier = np.einsum(einsums[n],e,k_s,field_fourier)*TWOPIIMG

  return np.fft.irfftn(curl_fourier,axes=(0,1,2),s=shapeFFT).reshape([N,n])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog option(s) [ASCIItable(s)]', description = """
Add column(s) containing curl of requested column(s).
Operates on periodic ordered three-dimensional data sets
of vector and tensor fields.
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
              3: {'name': 'vector',
                  'shape': [3],
                 },
              9: {'name': 'tensor',
                  'shape': [3,3],
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
    table.labels_append(['{}_curlFFT({})'.format(i+1,data['label']) 
                        for i in range(np.prod(np.array(data['shape'])))])                         # extend ASCII header with new labels
  table.head_write()

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray()
  grid,size = damask.util.coordGridAndSize(table.data[:,table.label_indexrange(options.pos)])

# ------------------------------------------ process value field -----------------------------------

  stack = [table.data]
  for data in active:
    # we need to reverse order here, because x is fastest,ie rightmost, but leftmost in our x,y,z notation
    stack.append(curlFFT(size[::-1],
                         table.data[:,table.label_indexrange(data['label'])].
                         reshape(grid[::-1].tolist()+data['shape'])))

# ------------------------------------------ output result -----------------------------------------

  if len(stack) > 1: table.data = np.hstack(tuple(stack))
  table.data_writeArray('%.12g')

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)
