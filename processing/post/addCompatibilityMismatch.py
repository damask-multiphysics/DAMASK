#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

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
          (idx // res[0]) % res[1], \
          (idx // res[0] // res[1]) % res[2] )

def index(location,res):

  return ( location[0] % res[0]                    + \
          (location[1] % res[1]) * res[0]          + \
          (location[2] % res[2]) * res[0] * res[1]   )        
# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options file[s]', description = """
Add column containing debug information
Operates on periodic ordered three-dimensional data sets.

""" + string.replace('$Id$','\n','\\n')
)


parser.add_option('--no-shape','-s',    dest='noShape', action='store_false', \
                                        help='do not calcuate shape mismatch [%default]')
parser.add_option('--no-volume','-v',   dest='noVolume', action='store_false', \
                                        help='do not calculate volume mismatch [%default]')
parser.add_option('-c','--coordinates', dest='coords', type='string',\
                                        help='column heading for coordinates [%default]')
parser.add_option('-f','--deformation', dest='F', action='extend', type='string', \
                                        help='heading(s) of columns containing deformation tensor values %default')

parser.set_defaults(noVolume = False)
parser.set_defaults(noShape = False)
parser.set_defaults(coords  = 'ip')
parser.set_defaults(F = 'f')

(options,filenames) = parser.parse_args()


datainfo = {                                                               # list of requested labels per datatype
             'F':     {'len':9,
                             'label':[]},
           }

if options.F != None:   datainfo['F']['label'] += options.F

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
  res = numpy.array([len(grid[0]),\
                     len(grid[1]),\
                     len(grid[2]),],'i')                                    # resolution is number of distinct coordinates found
  geomdim = res/numpy.maximum(numpy.ones(3,'d'),res-1.0)* \
              numpy.array([max(map(float,grid[0].keys()))-min(map(float,grid[0].keys())),\
                           max(map(float,grid[1].keys()))-min(map(float,grid[1].keys())),\
                           max(map(float,grid[2].keys()))-min(map(float,grid[2].keys())),\
                          ],'d')                                            # dimension from bounding box, corrected for cell-centeredness
  if res[2] == 1:
    geomdim[2] = min(geomdim[:2]/res[:2])

  N = res.prod()
  print '\t%s @ %s'%(geomdim,res)


# --------------- figure out columns to process

  key = '1_%s' %options.F
  if key not in table.labels:
    sys.stderr.write('column %s not found...\n'%key)
  else:
    F = numpy.array([0.0 for i in xrange(N*9)]).reshape([3,3]+list(res))
    if not options.noShape:  table.labels_append(['shapeMismatch(%s)' %options.F])
    if not options.noVolume: table.labels_append(['volMismatch(%s)'%options.F])
    column = table.labels.index(key)

# ------------------------------------------ assemble header ---------------------------------------  

  table.head_write()

# ------------------------------------------ read deformation gradient field -----------------------

  table.data_rewind()

  idx = 0
  while table.data_read():                                                  # read next data line of ASCII table
    (x,y,z) = location(idx,res)                                             # figure out (x,y,z) position from line count
    idx += 1
    F[0:3,0:3,x,y,z] = numpy.array(map(float,table.data[column:column+9]),'d').reshape(3,3)                                               
  
  Favg = damask.core.math.tensorAvg(F)

  if (res[0]%2 != 0 or res[1]%2 != 0 or (res[2] != 1 and res[2]%2 !=0)):
    print 'using linear reconstruction for uneven resolution'
    centres = damask.core.mesh.deformedCoordsLin(geomdim,F,Favg)
  else:
    centres = damask.core.mesh.deformedCoordsFFT(geomdim,F,1.0,Favg)
  
  nodes   = damask.core.mesh.nodesAroundCentres(geomdim,Favg,centres)
  if not options.noShape:   shapeMismatch = damask.core.mesh.shapeMismatch( geomdim,F,nodes,centres)
  if not options.noVolume: volumeMismatch = damask.core.mesh.volumeMismatch(geomdim,F,nodes)

# ------------------------------------------ process data ---------------------------------------  

  table.data_rewind()
  idx = 0
  while table.data_read():                                                  # read next data line of ASCII table
    (x,y,z) = location(idx,res)                                             # figure out (x,y,z) position from line count
    idx += 1
    if not options.noShape:  table.data_append( shapeMismatch[x,y,z])
    if not options.noVolume: table.data_append(volumeMismatch[x,y,z])
    
    table.data_write()                                                      # output processed line 

# ------------------------------------------ output result ---------------------------------------  

  table.output_flush()                                                      # just in case of buffered ASCII table

  file['input'].close()                                                     # close input ASCII table
  if file['name'] != 'STDIN':
    file['output'].close                                                    # close output ASCII table
    os.rename(file['name']+'_tmp',file['name'])                             # overwrite old one with tmp new
