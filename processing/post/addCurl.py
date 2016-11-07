#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

def curlFFT(geomdim,field):
 shapeFFT    = np.array(np.shape(field))[0:3]
 grid = np.array(np.shape(field)[2::-1])
 N = grid.prod()                                                                          # field size
 n = np.array(np.shape(field)[3:]).prod()                                                 # data size

 if   n == 3:   dataType = 'vector'
 elif n == 9:   dataType = 'tensor'

 field_fourier = np.fft.rfftn(field,axes=(0,1,2),s=shapeFFT)
 curl_fourier  = np.empty(field_fourier.shape,'c16')

# differentiation in Fourier space
 TWOPIIMG = 2.0j*math.pi
 k_sk = np.where(np.arange(grid[2])>grid[2]//2,np.arange(grid[2])-grid[2],np.arange(grid[2]))/geomdim[0]
 if grid[2]%2 == 0: k_sk[grid[2]//2] = 0                                                  # for even grid, set Nyquist freq to 0 (Johnson, MIT, 2011)
 
 k_sj = np.where(np.arange(grid[1])>grid[1]//2,np.arange(grid[1])-grid[1],np.arange(grid[1]))/geomdim[1]
 if grid[1]%2 == 0: k_sj[grid[1]//2] = 0                                                  # for even grid, set Nyquist freq to 0 (Johnson, MIT, 2011)
 
 k_si = np.arange(grid[0]//2+1)/geomdim[2]
 
 kk, kj, ki = np.meshgrid(k_sk,k_sj,k_si,indexing = 'ij')
 k_s = np.concatenate((ki[:,:,:,None],kj[:,:,:,None],kk[:,:,:,None]),axis = 3).astype('c16')
 
 e = np.zeros((3, 3, 3))
 e[0, 1, 2] = e[1, 2, 0] = e[2, 0, 1] = 1.0                                               # Levi-Civita symbols 
 e[0, 2, 1] = e[2, 1, 0] = e[1, 0, 2] = -1.0
 
 if dataType == 'tensor':                                                                 # tensor, 3x3 -> 3x3 
   curl_fourier = np.einsum('slm,ijkl,ijknm->ijksn',e,k_s,field_fourier)*TWOPIIMG
 elif dataType == 'vector':                                                               # vector, 3 -> 3
   curl_fourier = np.einsum('slm,ijkl,ijkm->ijks',e,k_s,field_fourier)*TWOPIIMG

 return np.fft.irfftn(curl_fourier,axes=(0,1,2),s=shapeFFT).reshape([N,n])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog option(s) [ASCIItable(s)]', description = """
Add column(s) containing curl of requested column(s).
Operates on periodic ordered three-dimensional data sets.
Deals with both vector- and tensor fields.

""", version = scriptID)

parser.add_option('-p','--pos','--periodiccellcenter',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'label of coordinates [%default]')
parser.add_option('-v','--vector',
                  dest = 'vector',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'label(s) of vector field values')
parser.add_option('-t','--tensor',
                  dest = 'tensor',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'label(s) of tensor field values')

parser.set_defaults(pos = 'pos',
                   )

(options,filenames) = parser.parse_args()

if options.vector is None and options.tensor is None:
  parser.error('no data column specified.')

# --- loop over input files ------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  items = {
            'tensor': {'dim': 9, 'shape': [3,3], 'labels':options.tensor, 'active':[], 'column': []},
            'vector': {'dim': 3, 'shape': [3],   'labels':options.vector, 'active':[], 'column': []},
          }
  errors  = []
  remarks = []
  column = {}
  
  if table.label_dimension(options.pos) != 3: errors.append('coordinates {} are not a vector.'.format(options.pos))
  else: colCoord = table.label_index(options.pos)

  for type, data in items.iteritems():
    for what in (data['labels'] if data['labels'] is not None else []):
      dim = table.label_dimension(what)
      if dim != data['dim']: remarks.append('column {} is not a {}.'.format(what,type))
      else:
        items[type]['active'].append(what)
        items[type]['column'].append(table.label_index(what))

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  for type, data in items.iteritems():
    for label in data['active']:
      table.labels_append(['{}_curlFFT({})'.format(i+1,label) for i in range(data['dim'])])         # extend ASCII header with new labels
  table.head_write()

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray()

  coords = [np.unique(table.data[:,colCoord+i]) for i in range(3)]
  mincorner = np.array(map(min,coords))
  maxcorner = np.array(map(max,coords))
  grid   = np.array(map(len,coords),'i')
  size   = grid/np.maximum(np.ones(3,'d'), grid-1.0) * (maxcorner-mincorner)                        # size from edge to edge = dim * n/(n-1) 
  size   = np.where(grid > 1, size, min(size[grid > 1]/grid[grid > 1]))                             # spacing for grid==1 equal to smallest among other ones

# ------------------------------------------ process value field -----------------------------------

  stack = [table.data]
  for type, data in items.iteritems():
    for i,label in enumerate(data['active']):
      # we need to reverse order here, because x is fastest,ie rightmost, but leftmost in our x,y,z notation
      stack.append(curlFFT(size[::-1],
                           table.data[:,data['column'][i]:data['column'][i]+data['dim']].
                           reshape(grid[::-1].tolist()+data['shape'])))

# ------------------------------------------ output result -----------------------------------------

  if len(stack) > 1: table.data = np.hstack(tuple(stack))
  table.data_writeArray('%.12g')

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)
