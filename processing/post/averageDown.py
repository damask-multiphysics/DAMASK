#!/usr/bin/env python

import os,re,sys,math,string,numpy,damask,time
from optparse import OptionParser, Option

# -----------------------------
class extendableOption(Option):
# -----------------------------
# used for definition of new option parser action 'extend', which enables to take multiple option arguments
# taken from online tutorial http://docs.python.org/library/optparse.html
  
  ACTIONS = Option.ACTIONS + ("extend",)
  STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
  TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
  ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

  def take_action(self, action, dest, opt, value, values, parser):
    if action == "extend":
      lvalue = value.split(",")
      values.ensure_value(dest, []).extend(lvalue)
    else:
      Option.take_action(self, action, dest, opt, value, values, parser)



def location(idx,res):

  return numpy.array([ idx  % res[0], \
                      (idx // res[0]) % res[1], \
                      (idx // res[0] // res[1]) % res[2] ])

def index(location,res):

  return ( location[0] % res[0]                    + \
          (location[1] % res[1]) * res[0]          + \
          (location[2] % res[2]) * res[0] * res[1]   )


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog [options] [file[s]]', description = """
Average each data block of size 'packing' into single values thus reducing the former resolution
to resolution/packing. (Requires numpy.)

""" + string.replace('$Id$','\n','\\n')
)

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
parser.set_defaults(coords     = 'ip')
parser.set_defaults(packing    = [2,2,2])
parser.set_defaults(shift      = [0,0,0])
parser.set_defaults(resolution = [0,0,0])
parser.set_defaults(dimension  = [0.0,0.0,0.0])

(options,filenames) = parser.parse_args()

if len(options.packing) < 3:
  parser.error('packing needs three parameters...')
if len(options.shift) < 3:
  parser.error('shift needs three parameters...')

options.packing = numpy.array(options.packing)
options.shift = numpy.array(options.shift)

prefix = 'averagedDown%ix%ix%i_'%(options.packing[0],options.packing[1],options.packing[2])
if numpy.any(options.shift != 0):
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
  

  try:
    locationCol = []
    for i,direction in enumerate(['x','y','z']):
      locationCol.append(table.labels.index('%s.%s'%(options.coords,direction))) # columns containing location data
    elemCol = table.labels.index('elem')                                    # columns containing location data
  except ValueError:
    print 'no coordinate data or element data found...'
    continue

  if (any(options.resolution)==0 or any(options.dimension)==0.0):
    grid = [{},{},{}]
    while table.data_read():                                                  # read next data line of ASCII table
      for j in range(3):
        grid[j][str(table.data[locationCol[j]])] = True                        # remember coordinate along x,y,z
    resolution = numpy.array([len(grid[0]),\
                              len(grid[1]),\
                              len(grid[2]),],'i')                             # resolution is number of distinct coordinates found
    dimension = resolution/numpy.maximum(numpy.ones(3,'d'),resolution-1.0)* \
                numpy.array([max(map(float,grid[0].keys()))-min(map(float,grid[0].keys())),\
                             max(map(float,grid[1].keys()))-min(map(float,grid[1].keys())),\
                             max(map(float,grid[2].keys()))-min(map(float,grid[2].keys())),\
                            ],'d')                                            # dimension from bounding box, corrected for cell-centeredness
    origin = numpy.array([min(map(float,grid[0].keys())),\
                          min(map(float,grid[1].keys())),\
                          min(map(float,grid[2].keys())),\
                         ],'d') - 0.5 * dimension / resolution
  else:
    resolution = numpy.array(options.resolution,'i')
    dimension = numpy.array(options.dimension,'d')
    origin = numpy.zeros(3,'d')

  if resolution[2] == 1:
    options.packing[2] = 1
    options.shift[2]   = 0
    dimension[2]       = min(dimension[:2]/resolution[:2])                    # z spacing equal to smaller of x or y spacing
  
  packing   = numpy.array(options.packing,'i')
  shift     = numpy.array(options.shift,'i')
  downSized = numpy.maximum(numpy.ones(3,'i'),resolution//packing)
  outSize   = numpy.ceil(numpy.array(resolution,'d')/numpy.array(packing,'d'))
  
  print '\t%s @ %s --> %s'%(dimension,resolution,downSized)
  
# ------------------------------------------ assemble header ---------------------------------------  

  table.head_write()

# ------------------------------------------ process data ---------------------------------------  
  table.data_rewind()
  data = numpy.zeros(outSize.tolist()+[len(table.labels)])
  p = numpy.zeros(3,'i')
  
  for p[2] in xrange(resolution[2]):
    for p[1] in xrange(resolution[1]):
      for p[0] in xrange(resolution[0]):
        d = ((p-shift)%resolution)//packing
        table.data_read()
        data[d[0],d[1],d[2],:] += numpy.array(table.data_asFloat(),'d')                        # convert to numpy array
   
  data /= packing.prod()


  elementSize = dimension/resolution*packing
  posOffset = (shift+[0.5,0.5,0.5])*elementSize
  elem = 1
  for c in xrange(downSized[2]):
    for b in xrange(downSized[1]):
      for a in xrange(downSized[0]):
        for i,x in enumerate([a,b,c]):
          data[a,b,c,locationCol[i]] = posOffset[i] + x*elementSize[i] + origin[i]
        data[a,b,c,elemCol] = elem
        table.data = data[a,b,c,:].tolist()
        table.data_write()                                                  # output processed line
        elem += 1

  
# ------------------------------------------ output result ---------------------------------------  

  table.output_flush()                                                      # just in case of buffered ASCII table

# ------------------------------------------ close file handles ---------------------------------------  

for file in files:
  file['input'].close()                                                     # close input ASCII table
  if file['name'] != 'STDIN':
    file['output'].close()                                                  # close output ASCII table
