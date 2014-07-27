#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,string
import numpy as np
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID = '$Id$'
scriptName = scriptID.split()[1]


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing divergence of requested column(s).
Operates on periodic ordered three-dimensional data sets.
Deals with both vector- and tensor-valued fields.

""", version = string.replace(scriptID,'\n','\\n')
)

accuracyChoices = ['2','4','6','8']
parser.add_option('--fdm',              dest='accuracy', action='extend', type='string', metavar='<int LIST>',
                                        help='degree of central difference accuracy (%s)'%(','.join(accuracyChoices)))
parser.add_option('--fft',              dest='fft', action='store_true',
                                        help='calculate divergence in Fourier space')
parser.add_option('-c','--coordinates', dest='coords', action='store', type='string', metavar = 'string',
                                        help='column heading for coordinates [%default]')
parser.add_option('-v','--vector',      dest='vector', action='extend', type='string', metavar='<string LIST>',
                                        help='heading of columns containing vector field values')
parser.add_option('-t','--tensor',      dest='tensor', action='extend', type='string', metavar='<string LIST>',
                                        help='heading of columns containing tensor field values')
parser.set_defaults(coords = 'ip')
parser.set_defaults(accuracy = [])
parser.set_defaults(fft = False)
parser.set_defaults(vector = [])
parser.set_defaults(tensor = [])

(options,filenames) = parser.parse_args()

if len(options.vector) + len(options.tensor) == 0:
  parser.error('no data column specified...')
if not set(options.accuracy).issubset(set(accuracyChoices)):
  parser.error('accuracy must be chosen from %s...'%(', '.join(accuracyChoices)))

if options.fft: options.accuracy.append('FFT')
if not options.accuracy:
  parser.error('no accuracy selected')

datainfo = {                                                               # list of requested labels per datatype
             'vector':     {'len':3,
                            'label':[]},
             'tensor':     {'len':9,
                            'label':[]},
           }

if options.vector != None:    datainfo['vector']['label'] += options.vector
if options.tensor != None:    datainfo['tensor']['label'] += options.tensor

# ------------------------------------------ setup file handles ------------------------------------
files = []
for name in filenames:
  if os.path.exists(name):
    files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],True)                                      # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.info_append(string.replace(scriptID,'\n','\\n') + '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out dimension and resolution ----------------------------------------------
  try:
    locationCol = table.labels.index('%s.x'%options.coords)                                         # columns containing location data
  except ValueError:
    file['croak'].write('no coordinate data found...\n'%key)
    continue

  grid = [{},{},{}]
  while table.data_read():                                                                          # read next data line of ASCII table
    for j in xrange(3):
      grid[j][str(table.data[locationCol+j])] = True                                                # remember coordinate along x,y,z
  resolution = np.array([len(grid[0]),\
                         len(grid[1]),\
                         len(grid[2]),],'i')                                                        # resolution is number of distinct coordinates found
  dimension = resolution/np.maximum(np.ones(3,'d'),resolution-1.0)* \
              np.array([max(map(float,grid[0].keys()))-min(map(float,grid[0].keys())),\
                        max(map(float,grid[1].keys()))-min(map(float,grid[1].keys())),\
                        max(map(float,grid[2].keys()))-min(map(float,grid[2].keys())),\
                        ],'d')                                                                      # dimension from bounding box, corrected for cell-centeredness
  if resolution[2] == 1:
    dimension[2] = min(dimension[:2]/resolution[:2])
  N = resolution.prod()
  
# --------------- figure out columns to process  --------------------------------------------------
  active     = defaultdict(list)
  column     = defaultdict(dict)
  values     = defaultdict(dict)
  divergence = defaultdict(dict)

  for datatype,info in datainfo.items():
    for label in info['label']:
      key = {True :'1_%s',
             False:'%s'   }[info['len']>1]%label
      if key not in table.labels:
        file['croak'].write('column %s not found...\n'%key)
      else:
        active[datatype].append(label)
        column[datatype][label] = table.labels.index(key)                                           # remember columns of requested data
        values[datatype][label] = np.array([0.0 for i in xrange(N*datainfo[datatype]['len'])]).\
                                           reshape(list(resolution)+[datainfo[datatype]['len']//3,3])
        if label not in divergence[datatype]: divergence[datatype][label] = {}
        for accuracy in options.accuracy:
          divergence[datatype][label][accuracy] = np.array([0.0 for i in xrange(N*datainfo[datatype]['len']//3)]).\
                                                           reshape(list(resolution)+[datainfo[datatype]['len']//3])
        
# ------------------------------------------ assemble header ---------------------------------------  
  for datatype,labels in active.items():                                                            # loop over vector,tensor
    for label in labels:
      for accuracy in options.accuracy:
        table.labels_append({True: ['%i_div%s(%s)'%(i+1,accuracy,label) for i in xrange(3)],       # extend ASCII header with new labels
                             False:['div%s(%s)'%(accuracy,label)]} [datatype == 'tensor'])
  table.head_write()

# ------------------------------------------ read value field --------------------------------------
  table.data_rewind()
  idx = 0
  while table.data_read():                                                                          # read next data line of ASCII table
    (x,y,z) = damask.util.gridLocation(idx,resolution)                                              # figure out (x,y,z) position from line count
    idx += 1
    for datatype,labels in active.items():                                                          # loop over vector,tensor
      for label in labels:                                                                          # loop over all requested curls
        values[datatype][label][x,y,z] = np.array(
                map(float,table.data[column[datatype][label]:
                                     column[datatype][label]+datainfo[datatype]['len']]),'d') \
                                     .reshape(datainfo[datatype]['len']//3,3)

# ------------------------------------------ process value field -----------------------------------
  for datatype,labels in active.items():                                                            # loop over vector,tensor
    for label in labels:                                                                            # loop over all requested divergencies
      for accuracy in options.accuracy:
        if accuracy == 'FFT':
          divergence[datatype][label][accuracy] = damask.core.math.divergenceFFT(dimension,values[datatype][label])
        else:
          divergence[datatype][label][accuracy] = damask.core.math.divergenceFDM(dimension,eval(accuracy)//2-1,values[datatype][label])

# ------------------------------------------ process data ---------------------------------------
  table.data_rewind()
  idx = 0
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    (x,y,z) = damask.util.gridLocation(idx,resolution)                                              # figure out (x,y,z) position from line count
    idx += 1
    for datatype,labels in active.items():                                                          # loop over vector,tensor
      for label in labels:                                                                          # loop over all requested 
        for accuracy in options.accuracy:
          table.data_append(list(divergence[datatype][label][accuracy][x,y,z].reshape(datainfo[datatype]['len']//3)))

    outputAlive = table.data_write()                                                                 # output processed line

# ------------------------------------------ output result ---------------------------------------  
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  file['input'].close()                                                                             # close input ASCII table (works for stdin)
  file['output'].close()                                                                            # close output ASCII table (works for stdout)
  os.rename(file['name']+'_tmp',file['name'])                                                       # overwrite old one with tmp new
