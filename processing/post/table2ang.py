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


parser.add_option('--coords',           dest='coords', type='string', \
                                        help='label of coords in ASCII table')
parser.add_option('--eulerangles',      dest='eulerangles', type='string', \
                                        help='label of euler angles in ASCII table')
parser.add_option('--defgrad',          dest='defgrad', type='string', \
                                        help='label of deformation gradient in ASCII table')
parser.add_option('-n','--normal',      dest='normal', type='float', nargs=3, \
                                        help='normal of slices to visualize')
parser.add_option('-u','--up',          dest='up', type='float', nargs=3,
                                        help='up direction of slices to visualize')
parser.add_option('-r','--resolution',  dest='res', type='int', nargs=3,
                                        help='up direction of slices to visualize')
parser.add_option('-c','--center',      dest='center', type='float', nargs=3,
                                        help='center of ang file in cube, negative for center')
parser.set_defaults(coords = 'coords')
parser.set_defaults(eulerangles = 'eulerangles')
parser.set_defaults(defgrad = 'f')
parser.set_defaults(normal = ['0.0','0.0','1.0'])
parser.set_defaults(up = ['1.0','0.0','0.0'])
parser.set_defaults(center = ['-1.0','-1.0','-1.0'])
parser.set_defaults(res = ['16','16','1'])
(options,filenames) = parser.parse_args()

datainfo = {
             'vector':     {'len':3,
                            'label':[]},
             'tensor':     {'len':9,
                            'label':[]}
           }
           
datainfo['vector']['label'].append(options.coords)
datainfo['vector']['label'].append(options.eulerangles)
datainfo['tensor']['label'].append(options.defgrad)

print options.res[0]
print options.res[1]
print options.res[2]

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
        active[datatype].append(label)
        column[datatype][label] = table.labels.index(key)                   # remember columns of requested data
        values[datatype][label] = numpy.array([0.0 for i in xrange(N*datainfo[datatype]['len'])])
        
# ------------------------------------------ read value field ---------------------------------------  

  table.data_rewind()
  idx = 0
  while table.data_read():                                                  # read next data line of ASCII table
    for datatype,labels in active.items():                                  # loop over vector,tensor
      for label in labels:                                                  # loop over all requested curls
        begin = idx*datainfo[datatype]['len']
        end = begin + datainfo[datatype]['len']
        values[datatype][label][begin:end]= numpy.array(map(float,table.data[column[datatype][label]:
                                             column[datatype][label]+datainfo[datatype]['len']]),'d')
    idx+=1
  hexagonal  = True
  if hexagonal: 
    scale = math.sin(1.0/3.0*math.pi)
  else: 
    scale = 1.0

  res0 = int(float(options.res[0])/scale)

  print 'res0', res0
  print  'res 1', options.res[1]
 
  if hexagonal: 
    NpointsSlice = res0//2*(int(options.res[1])-1)+(res0-res0//2)*int(options.res[1])
  else:
    NpointsSlice = res0*int(options.res[1])

  print NpointsSlice
  z = numpy.array(options.normal,dtype='float')
  z = z/numpy.linalg.norm(z)
  x = numpy.array(options.up,dtype='float')
  x = x/numpy.linalg.norm(x)
  y = numpy.cross(z,x)
  x = numpy.cross(y,z)
  x = x/numpy.linalg.norm(x)

  Favg = damask.core.math.tensorAvg(values['tensor']['%s'%(options.defgrad)].reshape(resolution[0],resolution[1],resolution[2],3,3))
  mySlice = numpy.zeros(NpointsSlice*3)
  eulerangles = values['vector']['%s'%options.eulerangles].reshape([3,N],order='F')
  for i in xrange(int(options.res[2])):
    idx = 0
    shift = 0
    offset =  numpy.array([0.5,0.5,0.5],dtype='float')/[float(options.res[0]),float(options.res[1]),float(options.res[2])]*[dimension[0],dimension[1],dimension[2]]
    for j in xrange(res0):
      if hexagonal: 
        res1=int(options.res[1])-j%2
        myOffset = offset +float(j%2)* numpy.array([0.0,0.5,0.0],dtype='float')/[float(options.res[0]),float(options.res[1]),float(options.res[2])]*[dimension[0],dimension[1],dimension[2]]
      else:
        res1=int(options.res[1])
        myOffset = offset
      for k in xrange(res1):
        mySlice[idx*3:idx*3+3] = numpy.dot(numpy.array([x,y,z],dtype='float'),
            numpy.array([j,k,i],dtype='float'))/[float(res0),float(options.res[0]),float(options.res[2])]*[dimension[0],dimension[1],dimension[2]]\
            + myOffset
        idx+=1
    mySlice = mySlice.reshape([3,NpointsSlice],order='F')
    indices=damask.core.math.math_nearestNeighborSearch(3,Favg,numpy.array(
      dimension,dtype='float'),NpointsSlice,N,mySlice,values['vector']['%s'%options.coords].reshape([3,N],order='F'))/27
    fileOut=open(os.path.join(os.path.dirname(name),os.path.splitext(os.path.basename(name))[0]+'_%s.ang'%i),'w')
    for line in getHeader(res0,res1,1.0):
      fileOut.write(line + '\n')
  
    # write data
    for idx in xrange(NpointsSlice):
      fileOut.write(''.join(['%10.5f'%positiveRadians(angle) for angle in eulerangles[:,indices[idx]]])+
           ' %10.5f %10.5f'%(mySlice[1,idx],mySlice[0,idx])+
           ' 100.0 1.0 0 1 1.0\n')

    fileOut.close()

