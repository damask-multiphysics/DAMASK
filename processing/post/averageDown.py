#!/usr/bin/env python

import os,re,sys,math,string,numpy,damask
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

parser.set_defaults(coords = 'ip')
parser.set_defaults(packing = [2,2,2])
parser.set_defaults(shift = [0,0,0])

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
    locationCol = table.labels.index('%s.x'%options.coords)                 # columns containing location data
  except ValueError:
    print 'no coordinate data found...'
    continue

  grid = [{},{},{}]
  while table.data_read():                                                  # read next data line of ASCII table
    for j in xrange(3):
      grid[j][str(table.data[locationCol+j])] = True                        # remember coordinate along x,y,z
  resolution = numpy.array([len(grid[0]),\
                            len(grid[1]),\
                            len(grid[2]),],'i')                             # resolution is number of distinct coordinates found
  dimension = resolution/numpy.maximum(numpy.ones(3,'d'),resolution-1.0)* \
              numpy.array([max(map(float,grid[0].keys()))-min(map(float,grid[0].keys())),\
                           max(map(float,grid[1].keys()))-min(map(float,grid[1].keys())),\
                           max(map(float,grid[2].keys()))-min(map(float,grid[2].keys())),\
                          ],'d')                                            # dimension from bounding box, corrected for cell-centeredness
  if resolution[2] == 1:
    options.packing[2] = 1
    options.shift[2]   = 0
    dimension[2]       = min(dimension[:2]/resolution[:2])

  downSized = numpy.maximum(numpy.ones(3,'i'),resolution//options.packing)

  print '\t%s @ %s --> %s'%(dimension,resolution,downSized)
  
# ------------------------------------------ assemble header ---------------------------------------  

  table.head_write()

# ------------------------------------------ process data ---------------------------------------  

  table.data_rewind()

  averagedDown = numpy.zeros(downSized.tolist()+[len(table.labels)])

  for z in xrange(-options.shift[2],-options.shift[2]+resolution[2]):
    for y in xrange(-options.shift[1],-options.shift[1]+resolution[1]):
      for x in xrange(-options.shift[0],-options.shift[0]+resolution[0]):
        table.data_read()
        data = numpy.array(table.data_asFloat(),'d')                        # convert to numpy array
        me = numpy.array((x,y,z),'i')                                       # my location as array
        data[locationCol:locationCol+3] -= dimension*(me//resolution)       # shift coordinates if periodic image
        (a,b,c) = (me%resolution)//options.packing                          # bin to condense my location into
        averagedDown[a,b,c,:] += data                                       # store the (coord-updated) data there

  averagedDown /= options.packing.prod()                                    # normalize data by element count

  for c in xrange(downSized[2]):
    for b in xrange(downSized[1]):
      for a in xrange(downSized[0]):
        table.data = averagedDown[a,b,c,:].tolist()
        table.data_write()                                                  # output processed line

  
# ------------------------------------------ output result ---------------------------------------  

  table.output_flush()                                                      # just in case of buffered ASCII table

# ------------------------------------------ close file handles ---------------------------------------  

for file in files:
  file['input'].close()                                                     # close input ASCII table
  if file['name'] != 'STDIN':
    file['output'].close()                                                  # close output ASCII table
