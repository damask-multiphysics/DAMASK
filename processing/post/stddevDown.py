#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
import damask
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [options] datafile[s]', description = """
Calculates the standard deviation of data in blocks of size 'packing' thus reducing the former resolution
to resolution/packing.

""", version = scriptID)

parser.add_option('-c','--coordinates', dest='coords', type='string',\
                                        help='column heading for coordinates [%default]')
parser.add_option('-p','--packing',     dest='packing', type='int', nargs=3, \
                                        help='dimension of packed group %default')
parser.add_option('-s','--shift',       dest='shift', type='int', nargs=3, \
                                        help='shift vector of packing stencil %default')
parser.add_option('-r','--resolution',  dest='resolution', type='int', nargs=3, \
                                        help='resolution in x,y,z [autodetect]')
parser.add_option('-d','--dimension',   dest='dimension', type='float', nargs=3, \
                                        help='dimension in x,y,z [autodetect]')
parser.set_defaults(coords     = 'ipinitialcoord')
parser.set_defaults(packing    = [2,2,2])
parser.set_defaults(shift      = [0,0,0])
parser.set_defaults(resolution = [0,0,0])
parser.set_defaults(dimension  = [0.0,0.0,0.0])

(options,filenames) = parser.parse_args()

if len(options.packing) < 3:
  parser.error('packing needs three parameters...')
if len(options.shift) < 3:
  parser.error('shift needs three parameters...')

options.packing = np.array(options.packing)
options.shift = np.array(options.shift)

prefix = 'stddevDown%ix%ix%i_'%(options.packing[0],options.packing[1],options.packing[2])
if np.any(options.shift != 0):
  prefix += 'shift%+i%+i%+i_'%(options.shift[0],options.shift[1],options.shift[2])

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    name = os.path.relpath(name)
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name),
                    'output':open(os.path.join(os.path.dirname(name),prefix+os.path.basename(name)),'w')})


# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': print file['name'],

  table = damask.ASCIItable(file['input'],file['output'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info
  table.info_append(string.replace('$Id$','\n','\\n') + \
                    '\t' + ' '.join(sys.argv[1:]))
  

# --------------- figure out size and grid ---------------------------------------------------------
  try:
    locationCol = table.labels.index('1_%s'%options.coords)                                         # columns containing location data
  except ValueError:
    try:
      locationCol = table.labels.index('%s.x'%options.coords)                                       # columns containing location data (legacy naming scheme)
    except ValueError:
      file['croak'].write('no coordinate data (1_%s/%s.x) found...\n'%(options.coords,options.coords))
      continue

  if (any(options.resolution)==0 or any(options.dimension)==0.0):
    grid = [{},{},{}]
    while table.data_read():                                                  # read next data line of ASCII table
      for j in xrange(3):
        grid[j][str(table.data[locationCol+j])] = True                        # remember coordinate along x,y,z
    resolution = np.array([len(grid[0]),\
                              len(grid[1]),\
                              len(grid[2]),],'i')                             # resolution is number of distinct coordinates found
    dimension = resolution/np.maximum(np.ones(3,'d'),resolution-1.0)* \
                np.array([max(map(float,grid[0].keys()))-min(map(float,grid[0].keys())),\
                             max(map(float,grid[1].keys()))-min(map(float,grid[1].keys())),\
                             max(map(float,grid[2].keys()))-min(map(float,grid[2].keys())),\
                            ],'d')                                            # dimension from bounding box, corrected for cell-centeredness
  else:
    resolution = np.array(options.resolution,'i')
    dimension = np.array(options.dimension,'d')

  if resolution[2] == 1:
    options.packing[2] = 1
    options.shift[2]   = 0
    dimension[2]       = min(dimension[:2]/resolution[:2])                    # z spacing equal to smaller of x or y spacing
  
  packing   = np.array(options.packing,'i')
  shift     = np.array(options.shift,'i')
  downSized = np.maximum(np.ones(3,'i'),resolution//packing)
  outSize   = np.ceil(np.array(resolution,'d')/np.array(packing,'d'))
  
  print '\t%s @ %s --> %s'%(dimension,resolution,downSized)
  
# ------------------------------------------ assemble header ---------------------------------------  

  table.head_write()

# ------------------------------------------ process data ---------------------------------------  

  dataavg = np.zeros(outSize.tolist()+[len(table.labels)])
  datavar = np.zeros(outSize.tolist()+[len(table.labels)])
  p = np.zeros(3,'i')
  
  table.data_rewind()
  for p[2] in xrange(resolution[2]):
    for p[1] in xrange(resolution[1]):
      for p[0] in xrange(resolution[0]):
        d = ((p-shift)%resolution)//packing
        table.data_read()
        dataavg[d[0],d[1],d[2],:] += np.array(table.data_asFloat(),'d')                        # convert to np array
   
  dataavg /= packing.prod()
  
  table.data_rewind()
  for p[2] in xrange(resolution[2]):
    for p[1] in xrange(resolution[1]):
      for p[0] in xrange(resolution[0]):
        d = ((p-shift)%resolution)//packing
        table.data_read()
        datavar[d[0],d[1],d[2],:] += (np.array(table.data_asFloat(),'d') - dataavg[d[0],d[1],d[2],:])**2  

  datavar = np.sqrt(datavar/packing.prod())
  
  posOffset = (shift+[0.5,0.5,0.5])*dimension/resolution
  elementSize = dimension/resolution*packing
  for c in xrange(downSized[2]):
    for b in xrange(downSized[1]):
      for a in xrange(downSized[0]):
        datavar[a,b,c,locationCol:locationCol+3] = posOffset + [a,b,c]*elementSize
        table.data = datavar[a,b,c,:].tolist()
        table.data_write()                                                  # output processed line
  
# ------------------------------------------ output result ---------------------------------------  

  table.output_flush()                                                      # just in case of buffered ASCII table

# ------------------------------------------ close file handles ---------------------------------------  

for file in files:
  table.input_close()                                                       # close input ASCII table
  if file['name'] != 'STDIN':
    table.output_close()                                                    # close output ASCII table
