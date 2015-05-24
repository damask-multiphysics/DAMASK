#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math,operator
import numpy as np
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

def divFFT(geomdim,field):
 grid = np.array(np.shape(field)[0:3])
 wgt = 1.0/np.array(grid).prod()

 field_fourier=np.fft.fftpack.rfftn(field,axes=(0,1,2))

 if len(np.shape(field)) == 4:
   dataType = 'vector'
   div_fourier=np.zeros(field_fourier.shape[0:3],'c16') # div is a scalar
 elif len(np.shape(field)) == 5:
   dataType = 'tensor'
   div_fourier=np.zeros(field_fourier.shape[0:4],'c16') # div is a vector

# differentiation in Fourier space
 k_s=np.zeros([3],'i')
 TWOPIIMG = (0.0+2.0j*math.pi)
 for i in xrange(grid[0]):
   k_s[0] = i
   if(i > grid[0]/2 ): k_s[0] = k_s[0] - grid[0]
   for j in xrange(grid[1]):
     k_s[1] = j
     if(j > grid[1]/2 ): k_s[1] = k_s[1] - grid[1]
     for k in xrange(grid[2]/2+1):
       k_s[2] = k
       if(k > grid[2]/2 ): k_s[2] = k_s[2] - grid[2]
       xi=np.array([k_s[2]/geomdim[2]+0.0j,k_s[1]/geomdim[1]+0.j,k_s[0]/geomdim[0]+0.j],'c16')
       if dataType == 'tensor': 
         for l in xrange(3):
           div_fourier[i,j,k,l] = sum(field_fourier[i,j,k,l,0:3]*xi) *TWOPIIMG
       elif dataType == 'vector': 
         div_fourier[i,j,k] = sum(field_fourier[i,j,k,0:3]*xi) *TWOPIIMG

 div=np.fft.fftpack.irfftn(div_fourier,axes=(0,1,2))
 print div.shape
 if dataType == 'tensor': 
   return  div.reshape([grid.prod(),3])
 if dataType == 'vector': 
   return  div.reshape([grid.prod(),1])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing divergence of requested column(s).
Operates on periodic ordered three-dimensional data sets.
Deals with both vector- and tensor-valued fields.

""", version = scriptID)

parser.add_option('-c','--coordinates', dest='coords', metavar='string',
                  help='column heading for coordinates [%default]')
parser.add_option('-v','--vector',      dest='vector', action='extend', metavar='<string LIST>',
                  help='heading of columns containing vector field values')
parser.add_option('-t','--tensor',      dest='tensor', action='extend', metavar='<string LIST>',
                  help='heading of columns containing tensor field values')
parser.set_defaults(coords = 'ipinitialcoord')

(options,filenames) = parser.parse_args()

if options.vector == None and options.tensor == None:
  parser.error('no data column specified...')

datainfo = {                                                                                         # list of requested labels per datatype
             'vector':     {'shape':[3],
                            'len':3,
                            'label':[]},
             'tensor':     {'shape':[3,3],
                            'len':9,
                            'label':[]},
           }

if options.vector != None:    datainfo['vector']['label'] = options.vector
if options.tensor != None:    datainfo['tensor']['label'] = options.tensor

# ------------------------------------------ setup file handles ------------------------------------
files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

#--- loop over input files -------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  table = damask.ASCIItable(file['input'],file['output'],False)                                     # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.data_readArray()

# --------------- figure out name of coordinate data (support for legacy .x notation) -------------
  coordLabels=['%i_%s'%(i+1,options.coords) for i in xrange(3)]                                     # store labels for column keys
  if not set(coordLabels).issubset(table.labels):
    directions = ['x','y','z']
    coordLabels=['%s.%s'%(options.coords,directions[i]) for i in xrange(3)]                         # store labels for column keys
    if not set(coordLabels).issubset(table.labels):
      file['croak'].write('no coordinate data (1_%s) found...\n'%options.coords)
      continue
  coordColumns = [table.labels.index(label) for label in coordLabels]

# --------------- figure out active columns -------------------------------------------------------
  active = defaultdict(list)
  for datatype,info in datainfo.items():
    for label in info['label']:
      key = '1_%s'%label
      if key not in table.labels:
        file['croak'].write('column %s not found...\n'%key)
      else:
        active[datatype].append(label)

# --------------- assemble new header (metadata and columns containing curl) ----------------------
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  for datatype,labels in active.items():                                                            # loop over vector,tensor
    for label in labels:
      table.labels_append(['divFFT(%s)'%(label) if datatype == 'vector' else
                           '%i_divFFT(%s)'%(i+1,label) for i in xrange(datainfo[datatype]['len']//3)])# extend ASCII header with new labels
  table.head_write()

# --------------- figure out size and grid ---------------------------------------------------------
  coords = [{},{},{}]
  for i in xrange(table.data.shape[0]):  
    for j in xrange(3):
      coords[j][str(table.data[i,coordColumns[j]])] = True
  grid = np.array(map(len,coords),'i')
  size = grid/np.maximum(np.ones(3,'d'),grid-1.0)* \
            np.array([max(map(float,coords[0].keys()))-min(map(float,coords[0].keys())),\
                      max(map(float,coords[1].keys()))-min(map(float,coords[1].keys())),\
                      max(map(float,coords[2].keys()))-min(map(float,coords[2].keys())),\
                      ],'d')                                                                        # size from bounding box, corrected for cell-centeredness
  for i, points in enumerate(grid):
    if points == 1:
      mask = np.ones(3,dtype=bool)
      mask[i]=0
      size[i] = min(size[mask]/grid[mask])                                                          # third spacing equal to smaller of other spacing 

# ------------------------------------------ process value field -----------------------------------
  div = defaultdict(dict)
  for datatype,labels in active.items():                                                            # loop over vector,tensor
    for label in labels:                                                                            # loop over all requested curls
      startColumn=table.labels.index('1_'+label)
      div[datatype][label]  = divFFT(size[::-1],                                                    # we need to reverse order here, because x is fastest,ie rightmost, but leftmost in our x,y,z notation
                              table.data[:,startColumn:startColumn+datainfo[datatype]['len']].\
                              reshape([grid[2],grid[1],grid[0]]+datainfo[datatype]['shape']))

# ------------------------------------------ add data ------------------------------------------
  for datatype,labels in active.items():                                                            # loop over vector,tensor
    for label in labels:                                                                            # loop over all requested curls
      for c in xrange(div[datatype][label][0,:].shape[0]):                                         # append column by column
        lastRow = table.data.shape[1]
        table.data=np.insert(table.data,lastRow,div[datatype][label][:,c],1)

# ------------------------------------------ output result -----------------------------------------
  table.data_writeArray('%.12g')
  table.input_close()                                                                               # close input ASCII table (works for stdin)
  table.output_close()                                                                              # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
