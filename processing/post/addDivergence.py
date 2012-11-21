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
  return ( idx  % res[0], \
         ( idx // res[0]) % res[1], \
         ( idx // res[0] // res[1]) % res[2] )

def index(location,res):
  return ( location[0] % res[0]                   + \
         ( location[1] % res[1]) * res[0]          + \
         ( location[2] % res[2]) * res[1] * res[0]   )



# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing divergence of requested column(s).
Operates on periodic ordered three-dimensional data sets.
Deals with both vector- and tensor-valued fields.

""" + string.replace('$Id$','\n','\\n')
)

accuracyChoices = ['2','4','6','8']

parser.add_option('--fdm',              dest='accuracy', action='extend', type='string', \
                                        help='degree of central difference accuracy (%s)'%(','.join(accuracyChoices)))
parser.add_option('--fft',              dest='fft', action='store_true', \
                                        help='calculate divergence in Fourier space')
parser.add_option('-c','--coordinates', dest='coords', type='string',\
                                        help='column heading for coordinates [%default]')
parser.add_option('-v','--vector',      dest='vector', action='extend', type='string', \
                                        help='heading of columns containing vector field values')
parser.add_option('-t','--tensor',      dest='tensor', action='extend', type='string', \
                                        help='heading of columns containing tensor field values')

parser.set_defaults(coords = 'ip')
parser.set_defaults(accuracy = [])
parser.set_defaults(fft = False)
parser.set_defaults(vector = [])
parser.set_defaults(tensor = [])

(options,filenames) = parser.parse_args()

if len(options.vector) + len(options.tensor) == 0:
  parser.error('no data column specified...')
  
for choice in options.accuracy:
  if choice not in accuracyChoices:
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

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w')})


# ------------------------------------------ loop over input files ---------------------------------------   

for file in files:
  if file['name'] != 'STDIN': print file['name'],

  table = damask.ASCIItable(file['input'],file['output'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info
  table.info_append(string.replace('$Id$','\n','\\n') + \
                    '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out dimension and resolution 
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
    dimension[2] = min(dimension[:2]/resolution[:2])

  N = resolution.prod()
  print '\t%s @ %s'%(dimension,resolution)

  
# --------------- figure out columns to process
  active = {}
  column = {}
  values = {}
  divergence   = {}

  head = []

  for datatype,info in datainfo.items():
    for label in info['label']:
      key = {True :'1_%s',
             False:'%s'   }[info['len']>1]%label
      if key not in table.labels:
        sys.stderr.write('column %s not found...\n'%key)
      else:
        if datatype not in active:     active[datatype] = []
        if datatype not in column:     column[datatype] = {}
        if datatype not in values:     values[datatype] = {}
        if datatype not in divergence: divergence[datatype] = {}
        if label not in divergence[datatype]: divergence[datatype][label] = {}
        active[datatype].append(label)
        column[datatype][label] = table.labels.index(key)                   # remember columns of requested data
        values[datatype][label] = numpy.array([0.0 for i in xrange(N*datainfo[datatype]['len'])]).\
                                           reshape(list(resolution)+[datainfo[datatype]['len']//3,3])
        for accuracy in options.accuracy:
          divergence[datatype][label][accuracy] = numpy.array([0.0 for i in xrange(N*datainfo[datatype]['len']//3)]).\
                                                           reshape(list(resolution)+[datainfo[datatype]['len']//3])
          if datatype == 'vector':                        # extend ASCII header with new labels
            table.labels_append(['div%s(%s)'%(accuracy,label)])
          if datatype == 'tensor':
            table.labels_append(['%i_div%s(%s)'%(i+1,accuracy,label) for i in xrange(3)])

        
# ------------------------------------------ assemble header ---------------------------------------  

  table.head_write()

# ------------------------------------------ read value field ---------------------------------------  

  table.data_rewind()

  idx = 0
  while table.data_read():                                                  # read next data line of ASCII table
    (x,y,z) = location(idx,resolution)                                     # figure out (x,y,z) position from line count
    idx += 1
    for datatype,labels in active.items():                                  # loop over vector,tensor
      for label in labels:                                                  # loop over all requested curls
        values[datatype][label][x,y,z] = numpy.array(
                                           map(float,table.data[column[datatype][label]:
                                                                column[datatype][label]+datainfo[datatype]['len']]),'d').reshape(datainfo[datatype]['len']//3,3)
# ------------------------------------------ process value field ---------------------------------------  

  for datatype,labels in active.items():                                  # loop over vector,tensor
    for label in labels:                                                  # loop over all requested divergencies
      for accuracy in options.accuracy:
        if accuracy == 'FFT':
          divergence[datatype][label][accuracy] = damask.core.math.divergenceFFT(dimension,values[datatype][label])
        else:
          divergence[datatype][label][accuracy] = damask.core.math.divergenceFDM(dimension,eval(accuracy)//2-1,values[datatype][label])
# ------------------------------------------ process data ---------------------------------------  

  table.data_rewind()
  idx = 0
  while table.data_read():                                                  # read next data line of ASCII table
    (x,y,z) = location(idx,resolution)                                     # figure out (x,y,z) position from line count
    idx += 1

    for datatype,labels in active.items():                                  # loop over vector,tensor
      for label in labels:                                                  # loop over all requested 
        for accuracy in options.accuracy:
          table.data_append(list(divergence[datatype][label][accuracy][x,y,z].reshape(datainfo[datatype]['len']//3)))

    table.data_write()                                                      # output processed line

  
# ------------------------------------------ output result ---------------------------------------  

  table.output_flush()                                                      # just in case of buffered ASCII table

  file['input'].close()                                                     # close input ASCII table
  if file['name'] != 'STDIN':
    file['output'].close                                                    # close output ASCII table
    os.rename(file['name']+'_tmp',file['name'])                             # overwrite old one with tmp new
