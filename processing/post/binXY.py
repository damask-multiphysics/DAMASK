#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = scriptID.split()[1][:-3]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Produces a binned grid of two columns from an ASCIItable, i.e. a two-dimensional probability density map.

""", version = scriptID)

parser.add_option('-d','--data',    dest='data', nargs=2, type='int', metavar='int int',
                                    help='columns containing x and y %default')
parser.add_option('-w','--weight',  dest='weight', metavar='int', type='int',
                                    help='column containing weight of (x,y) point [%default]')
parser.add_option('-b','--bins',    dest='bins', nargs=2, type='int', metavar='int int',
                                    help='number of bins in x and y direction %default')
parser.add_option('-t','--type',    dest='type', nargs=3, metavar='string string string',
                                    help='type (linear/log) of x, y, and z axis [linear]')
parser.add_option('-x','--xrange',  dest='xrange', nargs=2, type='float', metavar='float float',
                                    help='value range in x direction [auto]')
parser.add_option('-y','--yrange',  dest='yrange', nargs=2, type='float', metavar='float float',
                                    help='value range in y direction [auto]')
parser.add_option('-z','--zrange',  dest='zrange', nargs=2, type='float', metavar='float float',
                                    help='value range in z direction [auto]')
parser.add_option('-i','--invert',  dest='invert', action='store_true',
                                    help='invert probability density [%default]')

parser.set_defaults(data = (1,2))
parser.set_defaults(weight = None)
parser.set_defaults(bins = (10,10))
parser.set_defaults(type = ('linear','linear','linear'))
parser.set_defaults(xrange = (0.0,0.0))
parser.set_defaults(yrange = (0.0,0.0))
parser.set_defaults(zrange = (0.0,0.0))
parser.set_defaults(invert = False)

(options,filenames) = parser.parse_args()

range =  np.array([np.array(options.xrange),
                   np.array(options.yrange),
                   np.array(options.zrange)])
grid =   np.zeros(options.bins,'i')
result = np.zeros((options.bins[0]*options.bins[1],3),'f')

prefix='binned%i-%i_'%(options.data[0],options.data[1])+ \
                      ('weighted%i_'%(options.weight) if options.weight != None else '')

# ------------------------------------------ setup file handles ------------------------------------
files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

# ------------------------------------------ loop over input files ---------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  skip = int(file['input'].readline().split()[0])
  for i in xrange(skip): headers = file['input'].readline().split()
  data = np.loadtxt(file['input'],usecols=np.array(options.data+((options.weight,) if options.weight != None else ()))-1)
  file['input'].close()                                                                              # close input ASCII table

  for i in (0,1):                                                                                    # check data range for x and y
    if (range[i] == 0.0).all(): range[i] = [data[:,i].min(),data[:,i].max()]
    if options.type[i].lower() == 'log':                                                             # if log scale
      data[:,i] = np.log(data[:,i])                                                                  # change x,y coordinates to log
      range[i] = np.log(range[i])                                                                    # change range to log, too
  
  delta = range[:,1]-range[:,0]
  
  for i in xrange(len(data)):
    x = int(options.bins[0]*(data[i,0]-range[0,0])/delta[0])
    y = int(options.bins[1]*(data[i,1]-range[1,0])/delta[1])
    if x >=0 and x < options.bins[0] and y >= 0 and y < options.bins[1]: grid[x,y] += 1 if options.weight == None else data[i,2]
  
  if (range[2] == 0.0).all(): range[2] = [grid.min(),grid.max()]
  if (range[2] == 0.0).all():                                                                       # no data in grid?
    file['croak'].write('no data found on grid...\n')
    range[2,:] = np.array([0.0,1.0])                                                                # making up arbitrary z range
  if options.type[2].lower() == 'log':
    grid = np.log(grid)
    range[2] = np.log(range[2])
    
  delta[2] = range[2,1]-range[2,0]

  
  i = 0
  for x in xrange(options.bins[0]):
    for y in xrange(options.bins[1]):
      result[i,:] = [range[0,0]+delta[0]/options.bins[0]*(x+0.5),
                     range[1,0]+delta[1]/options.bins[1]*(y+0.5),
                     min(1.0,max(0.0,(grid[x,y]-range[2,0])/delta[2]))]
      if options.type[0].lower() == 'log': result[i,0] = np.exp(result[i,0])
      if options.type[1].lower() == 'log': result[i,1] = np.exp(result[i,1])
      if options.invert:                   result[i,2] = 1.0-result[i,2]
      i += 1

 # ------------------------------------------ output result -----------------------------------------       
  file['output'].write('1\thead\n')
  file['output'].write('bin_%s\tbin_%s\tz\n'%(headers[options.data[0]-1],headers[options.data[1]-1]))
  np.savetxt(file['output'],result)
  file['output'].close()                                                                            # close output ASCII table
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',\
              os.path.join(os.path.dirname(file['name']),prefix+os.path.basename(file['name'])))
