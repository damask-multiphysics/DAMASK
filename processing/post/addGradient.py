#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

def gradFFT(geomdim,field):
 grid = np.array(np.shape(field)[2::-1])
 N = grid.prod()                                                                                    # field size
 n = np.array(np.shape(field)[3:]).prod()                                                           # data size
 if   n == 3:   dataType = 'vector'
 elif n == 1:   dataType = 'scalar'

 field_fourier = np.fft.fftpack.rfftn(field,axes=(0,1,2))
 grad_fourier   = np.zeros(field_fourier.shape+(3,),'c16')

# differentiation in Fourier space
 k_s = np.zeros([3],'i')
 TWOPIIMG = 2.0j*math.pi
 for i in xrange(grid[2]):
   k_s[0] = i
   if grid[2]%2 == 0 and i == grid[2]//2:  k_s[0] = 0                                                            # for even grid, set Nyquist freq to 0 (Johnson, MIT, 2011)
   elif i > grid[2]//2:                    k_s[0] -= grid[2]

   for j in xrange(grid[1]):
     k_s[1] = j
     if grid[1]%2 == 0 and j == grid[1]//2: k_s[1] = 0                                                          # for even grid, set Nyquist freq to 0 (Johnson, MIT, 2011)
     elif j > grid[1]//2:                   k_s[1] -= grid[1]

     for k in xrange(grid[0]//2+1):
       k_s[2] = k
       if grid[0]%2 == 0 and k == grid[0]//2: k_s[2] = 0                                                       # for even grid, set Nyquist freq to 0 (Johnson, MIT, 2011)

       xi = (k_s/geomdim)[2::-1].astype('c16')                                                                 # reversing the field order
       
       grad_fourier[i,j,k,0,:] = field_fourier[i,j,k,0]*xi *TWOPIIMG                                           # vector field from scalar data

       if dataType == 'vector':
         grad_fourier[i,j,k,1,:] = field_fourier[i,j,k,1]*xi *TWOPIIMG                                         # tensor field from vector data
         grad_fourier[i,j,k,2,:] = field_fourier[i,j,k,2]*xi *TWOPIIMG

 return np.fft.fftpack.irfftn(grad_fourier,axes=(0,1,2)).reshape([N,3*n])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing gradient of requested column(s).
Operates on periodic ordered three-dimensional data sets.
Deals with both vector- and scalar fields.

""", version = scriptID)

parser.add_option('-c','--coordinates',
                  dest = 'coords',
                  type = 'string', metavar='string',
                  help = 'column heading for coordinates [%default]')
parser.add_option('-v','--vector',
                  dest = 'vector',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'heading of columns containing vector field values')
parser.add_option('-s','--scalar',
                  dest = 'scalar',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'heading of columns containing scalar field values')

parser.set_defaults(coords = 'ipinitialcoord',
                   )

(options,filenames) = parser.parse_args()

if options.vector is None and options.scalar is None:
  parser.error('no data column specified.')

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  items = {
            'scalar': {'dim': 1, 'shape': [1], 'labels':options.scalar, 'active':[], 'column': []},
            'vector': {'dim': 3, 'shape': [3],   'labels':options.vector, 'active':[], 'column': []},
          }
  errors  = []
  remarks = []
  column = {}
  
  if table.label_dimension(options.coords) != 3: errors.append('coordinates {} are not a vector.'.format(options.coords))
  else: colCoord = table.label_index(options.coords)

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
      table.labels_append(['{}_gradFFT({})'.format(i+1,label) for i in xrange(3 * data['dim'])])        # extend ASCII header with new labels # grad increases the field dimension by one
  table.head_write()

# --------------- figure out size and grid ---------------------------------------------------------

  table.data_readArray()

  coords = [np.unique(table.data[:,colCoord+i]) for i in xrange(3)]
  mincorner = np.array(map(min,coords))
  maxcorner = np.array(map(max,coords))
  grid   = np.array(map(len,coords),'i')
  size   = grid/np.maximum(np.ones(3,'d'), grid-1.0) * (maxcorner-mincorner)                        # size from edge to edge = dim * n/(n-1) 
  size   = np.where(grid > 1, size, min(size[grid > 1]/grid[grid > 1]))     

# ------------------------------------------ process value field -----------------------------------

  stack = [table.data]
  for type, data in items.iteritems():
    for i,label in enumerate(data['active']):
      stack.append(gradFFT(size[::-1],                                                              # we need to reverse order here, because x is fastest,ie rightmost, but leftmost in our x,y,z notation
                           table.data[:,data['column'][i]:data['column'][i]+data['dim']].
                           reshape([grid[2],grid[1],grid[0]]+data['shape'])))

# ------------------------------------------ output result -----------------------------------------

  if len(stack) > 1: table.data = np.hstack(tuple(stack))
  table.data_writeArray('%.12g')

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)
