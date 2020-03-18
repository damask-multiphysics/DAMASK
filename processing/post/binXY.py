#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import numpy as np

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Produces a binned grid of two columns from an ASCIItable, i.e. a two-dimensional probability density map.

""", version = scriptID)

parser.add_option('-d','--data',
                  dest = 'data',
                  type = 'string', nargs = 2, metavar = 'string string',
                  help = 'column labels containing x and y ')
parser.add_option('-w','--weight',
                  dest = 'weight',
                  type = 'string', metavar = 'string',
                  help = 'column label containing weight of (x,y) point')
parser.add_option('-b','--bins',
                  dest = 'bins',
                  type = 'int', nargs = 2, metavar = 'int int',
                  help = 'number of bins in x and y direction [%default]')
parser.add_option('-t','--type',
                  dest = 'type',
                  type = 'string', nargs = 3, metavar = 'string string string',
                  help = 'type (linear/log) of x, y, and z axis [%default]')
parser.add_option('-x','--xrange',
                  dest = 'xrange',
                  type = 'float', nargs = 2, metavar = 'float float',
                  help = 'min max limits in x direction (optional)')
parser.add_option('-y','--yrange',
                  dest = 'yrange',
                  type = 'float', nargs = 2, metavar = 'float float',
                  help = 'min max limits in y direction (optional)')
parser.add_option('-z','--zrange',
                  dest = 'zrange',
                  type = 'float', nargs = 2, metavar = 'float float',
                  help = 'min max limits in z direction (optional)')
parser.add_option('-i','--invert',
                  dest = 'invert',
                  action = 'store_true',
                  help = 'invert probability density')
parser.add_option('-r','--rownormalize',
                  dest = 'normRow',
                  action = 'store_true',
                  help = 'normalize probability density in each row')
parser.add_option('-c','--colnormalize',
                  dest = 'normCol',
                  action = 'store_true',
                  help = 'normalize probability density in each column')

parser.set_defaults(bins = (10,10),
                    type = ('linear','linear','linear'),
                    xrange = (0.0,0.0),
                    yrange = (0.0,0.0),
                    zrange = (0.0,0.0),
                   )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

minmax = np.array([options.xrange,options.yrange,options.zrange])
result = np.empty((options.bins[0],options.bins[1],3),'f')

if options.data is None: parser.error('no data columns specified.')

for name in filenames:
  damask.util.report(scriptName,name)

  table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
  data = np.hstack((table.get(options.data[0]),table.get(options.data[1])))

  for c in (0,1):                                                                                   # check data minmax for x and y (i = 0 and 1)
    if (minmax[c] == 0.0).all(): minmax[c] = [data[:,c].min(),data[:,c].max()]
    if options.type[c].lower() == 'log':                                                            # if log scale
      data[:,c] = np.log(data[:,c])                                                                 # change x,y coordinates to log
      minmax[c] = np.log(minmax[c])                                                                 # change minmax to log, too

  delta = minmax[:,1]-minmax[:,0]
  (grid,xedges,yedges) = np.histogram2d(data[:,0],data[:,1],
                                        bins=options.bins,
                                        range=minmax[:2],
                                        weights=table.get(options.weight) if options.weight else None)
  if options.normCol:
    for x in range(options.bins[0]):
      sum = np.sum(grid[x,:])
      if sum > 0.0:
        grid[x,:] /= sum
  if options.normRow:
    for y in range(options.bins[1]):
      sum = np.sum(grid[:,y])
      if sum > 0.0:
        grid[:,y] /= sum

  if (minmax[2] == 0.0).all(): minmax[2] = [grid.min(),grid.max()]                                   # auto scale from data
  if minmax[2,0] == minmax[2,1]:
    minmax[2,0] -= 1.
    minmax[2,1] += 1.
  if (minmax[2] == 0.0).all():                                                                       # no data in grid?
    damask.util.croak('no data found on grid...')
    minmax[2,:] = np.array([0.0,1.0])                                                                # making up arbitrary z minmax
  if options.type[2].lower() == 'log':
    grid = np.log(grid)
    minmax[2] = np.log(minmax[2])

  delta[2] = minmax[2,1]-minmax[2,0]

  for x in range(options.bins[0]):
    for y in range(options.bins[1]):
      result[x,y,:] = [minmax[0,0]+delta[0]/options.bins[0]*(x+0.5),
                       minmax[1,0]+delta[1]/options.bins[1]*(y+0.5),
                       np.clip((grid[x,y]-minmax[2,0])/delta[2],0.0,1.0)]

  for c in (0,1):
    if options.type[c].lower() == 'log': result[:,:,c] = np.exp(result[:,:,c])

  if options.invert: result[:,:,2] = 1.0 - result[:,:,2]

  comments = scriptID + '\t' + ' '.join(sys.argv[1:])
  shapes = {'bin_%s'%options.data[0]:(1,),'bin_%s'%options.data[1]:(1,),'z':(1,)}
  table  = damask.Table(result.reshape(options.bins[0]*options.bins[1],3),shapes,[comments])
  if name:
    outname = os.path.join(os.path.dirname(name),'binned-{}-{}_'.format(*options.data) +
                                                ('weighted-{}_'.format(options.weight) if options.weight else '') +
                                                 os.path.basename(name))
    table.to_ASCII(outname)
  else:
    table.to_ASCII(sys.stdout)
