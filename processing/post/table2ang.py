#!/usr/bin/env python

import os,re,sys,math,string,damask,numpy
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

def integerFactorization(i):
  
  j = int(math.floor(math.sqrt(float(i))))
  while (j>1 and int(i)%j != 0):
    j -= 1
  return j

def positiveRadians(angle):

  angle = math.radians(float(angle))
  while angle < 0.0:
    angle += 2.0*math.pi

  return angle


def getHeader(sizeX,sizeY,step):
  
  return [ \
  '# TEM_PIXperUM          1.000000', \
  '# x-star                0.509548', \
  '# y-star                0.795272', \
  '# z-star                0.611799', \
  '# WorkingDistance       18.000000', \
  '#', \
  '# Phase                 1', \
  '# MaterialName          Al', \
  '# Formula               Fe', \
  '# Info', \
  '# Symmetry              43', \
  '# LatticeConstants      2.870 2.870 2.870  90.000  90.000  90.000', \
  '# NumberFamilies        4', \
  '# hklFamilies           1  1  0 1 0.000000 1', \
  '# hklFamilies           2  0  0 1 0.000000 1', \
  '# hklFamilies           2  1  1 1 0.000000 1', \
  '# hklFamilies           3  1  0 1 0.000000 1', \
  '# Categories            0 0 0 0 0 ', \
  '#', \
  '# GRID: SquareGrid', \
  '# XSTEP: ' + str(step), \
  '# YSTEP: ' + str(step), \
  '# NCOLS_ODD: ' + str(sizeX), \
  '# NCOLS_EVEN: ' + str(sizeX), \
  '# NROWS: ' + str(sizeY), \
  '#', \
  '# OPERATOR: ODFsammpling', \
  '#', \
  '# SAMPLEID: ', \
  '#', \
  '# SCANID: ', \
  '#', \
  ]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Builds an ang file out of ASCII table.

""" + string.replace('$Id$','\n','\\n')
)


parser.add_option('--coords',            dest='coords', type='string', \
                                        help='label of coords in ASCII table')
parser.set_defaults(norm = 'ip')


(options,filenames) = parser.parse_args()

datainfo = {
             'vector':     {'len':3,
                            'label':[]}
           }
           
datainfo['vector']['label'] += ['eulerangles']
# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name)})


# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': print file['name']

  table = damask.ASCIItable(file['input'])                                  # open ASCII_table for reading
  table.head_read()                                                         # read ASCII header info

# --------------- figure out dimension and resolution 
  try:
    locationCol = table.labels.index('ip.x')                                # columns containing location data
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
        values[datatype][label] = numpy.array([0.0 for i in xrange(N*datainfo[datatype]['len'])])
        
# ------------------------------------------ read value field ---------------------------------------  

  table.data_rewind()
  print values['vector']['eulerangles']
  print numpy.shape(values['vector']['eulerangles'])
  idx = 0
  while table.data_read():                                                  # read next data line of ASCII table
    for datatype,labels in active.items():                                  # loop over vector,tensor
      for label in labels:                                                  # loop over all requested curls
        values[datatype][label][idx:idx+3]= numpy.array(map(float,table.data[column[datatype][label]:
                                             column[datatype][label]+datainfo[datatype]['len']]),'d')
        idx+=3
  for z in xrange(resolution[2]):
    fileOut=open(os.path.join(os.path.dirname(name),os.path.splitext(os.path.basename(name))[0]+'_%s.ang'%z),'w')
    for line in getHeader(resolution[0],resolution[1],1.0):
      fileOut.write(line + '\n')
    
    # write data
    for counter in xrange(resolution[0]*resolution[1]):
      print z*resolution[0]*resolution[1]*3+counter*3,' - ',z*resolution[0]*resolution[1]*3+counter*3+3
      fileOut.write(''.join(['%10.5f'%positiveRadians(angle) for angle in values['vector']['eulerangles'][z*resolution[0]*resolution[1]*3+counter*3:z*resolution[0]*resolution[1]*3+counter*3+3]])+
           ''.join(['%10.5f'%coord for coord in [counter%resolution[0],counter//resolution[1]]])+
           ' 100.0 1.0 0 1 1.0\n')
      counter += 1
    
    fileOut.close()

