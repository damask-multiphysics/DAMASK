#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Produces a binned grid of two columns from an ASCIItable, i.e. a two-dimensional probability density map.

""", version = scriptID)

parser.add_option('-d','--data',    dest='data', nargs=2, type='string', metavar='string string',
                                    help='column labels containing x and y %default')
parser.add_option('-w','--weight',  dest='weight', metavar='string', type='string',
                                    help='column label containing weight of (x,y) point [%default]')
parser.add_option('-b','--bins',    dest='bins', nargs=2, type='int', metavar='int int',
                                    help='number of bins in x and y direction %default')
parser.add_option('-t','--type',    dest='type', nargs=3, metavar='string string string',
                                    help='type (linear/log) of x, y, and z axis [linear]')
parser.add_option('-x','--xrange',  dest='xrange', nargs=2, type='float', metavar='float float',
                                    help='value minmax in x direction [auto]')
parser.add_option('-y','--yrange',  dest='yrange', nargs=2, type='float', metavar='float float',
                                    help='value minmax in y direction [auto]')
parser.add_option('-z','--zrange',  dest='zrange', nargs=2, type='float', metavar='float float',
                                    help='value minmax in z direction [auto]')
parser.add_option('-i','--invert',  dest='invert', action='store_true',
                                    help='invert probability density [%default]')
parser.add_option('-r','--rownormalize',  dest='normRow', action='store_true',
                                    help='normalize probability density in each row [%default]')
parser.add_option('-c','--colnormalize',  dest='normCol', action='store_true',
                                    help='normalize probability density in each column [%default]')

parser.set_defaults(data = None)
parser.set_defaults(weight = None)
parser.set_defaults(bins = (10,10))
parser.set_defaults(type = ('linear','linear','linear'))
parser.set_defaults(xrange = (0.0,0.0))
parser.set_defaults(yrange = (0.0,0.0))
parser.set_defaults(zrange = (0.0,0.0))
parser.set_defaults(invert = False)
parser.set_defaults(normRow = False)
parser.set_defaults(normCol = False)

(options,filenames) = parser.parse_args()

minmax =  np.array([np.array(options.xrange),
                    np.array(options.yrange),
                    np.array(options.zrange)])
grid =   np.zeros(options.bins,'f')
result = np.zeros((options.bins[0],options.bins[1],3),'f')

datainfo = {                                                                                        # list of requested labels per datatype
             'scalar':     {'len':1,
                            'label':[]},
           }

if options.data != None:   datainfo['scalar']['label'] += options.data
if options.weight != None: datainfo['scalar']['label'] += options.weight

if len(datainfo['scalar']['label']) < 2:
  parser.error('missing column labels')

# --- loop over input files -------------------------------------------------------------------------
if filenames == []:
  filenames = ['STDIN']

for name in filenames:
  if name == 'STDIN':
    file = {'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m\n')
  else:
    if not os.path.exists(name): continue
    file = {'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],buffered = False)                          # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info

# --------------- figure out columns to process  ---------------------------------------------------
  active = []
  column = {}

  for label in datainfo['scalar']['label']:
    if label in table.labels:
      active.append(label)
      column[label] = table.labels.index(label)                                                     # remember columns of requested data
    else:
      file['croak'].write('column %s not found...\n'%label)
       
# ------------------------------------------ assemble header ---------------------------------------
  table.info_clear()
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels = ['bin_%s'%options.data[0],'bin_%s'%options.data[1],'z']
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  table.data_readArray([column[label] for label in active])

  for i in (0,1):                                                                                    # check data minmax for x and y
    if (minmax[i] == 0.0).all(): minmax[i] = [table.data[:,i].min(),table.data[:,i].max()]
    if options.type[i].lower() == 'log':                                                             # if log scale
      table.data[:,i] = np.log(table.data[:,i])                                                                  # change x,y coordinates to log
      minmax[i] = np.log(minmax[i])                                                                    # change minmax to log, too

  delta = minmax[:,1]-minmax[:,0]
  
  for i in xrange(len(table.data)):
    x = int(options.bins[0]*(table.data[i,0]-minmax[0,0])/delta[0])
    y = int(options.bins[1]*(table.data[i,1]-minmax[1,0])/delta[1])
    if x >= 0 and x < options.bins[0] and y >= 0 and y < options.bins[1]:
      grid[x,y] += 1. if options.weight == None else table.data[i,2]                                        # count (weighted) occurrences

  if options.normCol:
    for x in xrange(options.bins[0]):
      sum = np.sum(grid[x,:])
      if sum > 0.0:
        grid[x,:] /= sum
  if options.normRow:
    for y in xrange(options.bins[1]):
      sum = np.sum(grid[:,y])
      if sum > 0.0:
        grid[:,y] /= sum
  
  if (minmax[2] == 0.0).all(): minmax[2] = [grid.min(),grid.max()]                                   # auto scale from data
  if minmax[2,0] == minmax[2,1]:
    minmax[2,0] -= 1.
    minmax[2,1] += 1.
  if (minmax[2] == 0.0).all():                                                                       # no data in grid?
    file['croak'].write('no data found on grid...\n')
    minmax[2,:] = np.array([0.0,1.0])                                                                # making up arbitrary z minmax
  if options.type[2].lower() == 'log':
    grid = np.log(grid)
    minmax[2] = np.log(minmax[2])
    
  delta[2] = minmax[2,1]-minmax[2,0]

  for x in xrange(options.bins[0]):
    for y in xrange(options.bins[1]):
      result[x,y,:] = [minmax[0,0]+delta[0]/options.bins[0]*(x+0.5),
                       minmax[1,0]+delta[1]/options.bins[1]*(y+0.5),
                       min(1.0,max(0.0,(grid[x,y]-minmax[2,0])/delta[2]))]

  for i in xrange(2):
    if options.type[i].lower() == 'log': result[:,:,i] = np.exp(result[:,:,i])

  if options.invert: result[:,:,2] = 1.0 - result[:,:,2]

 # ------------------------------------------ output result -----------------------------------------       
  prefix = 'binned%s-%s_'%(options.data[0],options.data[1])+ \
                          ('weighted%s_'%(options.weight) if options.weight != None else '')
  np.savetxt(file['output'],result.reshape(options.bins[0]*options.bins[1],3))
  file['output'].close()                                                                            # close output ASCII table
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',\
              os.path.join(os.path.dirname(file['name']),prefix+os.path.basename(file['name'])))
