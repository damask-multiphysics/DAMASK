#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,numpy
from optparse import OptionParser, Option

scriptID = '$Id$'
scriptName = scriptID.split()[1]

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



# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Produces a binned grid of two columns from an ASCIItable, i.e. a two-dimensional probability density map.
""" + string.replace(scriptID,'\n','\\n')
)


parser.add_option('-d','--data',    dest='data', nargs=2, type='int',
                                    help='columns containing x and y')
parser.add_option('-w','--weight',  dest='weight', type='int',
                                    help='column containing weight of (x,y) point')
parser.add_option('-b','--bins',    dest='bins', nargs=2, type='int',
                                    help='number of bins in x and y direction')
parser.add_option('-t','--type',    dest='type', nargs=3, type='string',
                                    help='type of x, y, and z axis [linear]')
parser.add_option('-x','--xrange',  dest='xrange', nargs=2, type='float',
                                    help='value range in x direction [auto]')
parser.add_option('-y','--yrange',  dest='yrange', nargs=2, type='float',
                                    help='value range in y direction [auto]')
parser.add_option('-z','--zrange',  dest='zrange', nargs=2, type='float',
                                    help='value range in z direction [auto]')
parser.add_option('-i','--invert',  dest='invert', action='store_true',
                                    help='invert probability density')

parser.set_defaults(data = [1,2])
parser.set_defaults(weight = None)
parser.set_defaults(bins = [10,10])
parser.set_defaults(type = ['linear','linear','linear'])
parser.set_defaults(xrange = [0.0,0.0])
parser.set_defaults(yrange = [0.0,0.0])
parser.set_defaults(zrange = [0.0,0.0])
parser.set_defaults(invert = False)

(options,filenames) = parser.parse_args()

range =  numpy.array([numpy.array(options.xrange),
                      numpy.array(options.yrange),
                      numpy.array(options.zrange)])
grid =   numpy.zeros(options.bins,'i')
result = numpy.zeros((options.bins[0]*options.bins[1],3),'f')
# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':   'STDIN',
                'input':  sys.stdin,
                'output': sys.stdout,
                'croak':  sys.stderr,
               })
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':   name,
                    'input':  open(name),
                    'output': open(os.path.splitext(name)[0]+ \
                                   '_binned%i-%i'%(options.data[0],options.data[1])+ \
                                   ('_weighted%i'%(options.weight) if options.weight != None else '')+ \
                                   os.path.splitext(name)[1],'w'),
                    'croak':  sys.stderr,
                    })

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  skip = int(file['input'].readline().split()[0])
  for i in xrange(skip): headers = file['input'].readline().split()
  data = numpy.loadtxt(file['input'],usecols=numpy.array(options.data+((options.weight,) if options.weight != None else ()))-1)
  file['input'].close()                                                   # close input ASCII table

  for i in (0,1):                                                         # check data range for x and y
    if (range[i] == 0.0).all(): range[i] = [data[:,i].min(),data[:,i].max()]
    if options.type[i].lower() == 'log':                                  # if log scale
      data[:,i] = numpy.log(data[:,i])                                    # change x,y coordinates to log
      range[i] = numpy.log(range[i])                                      # change range to log, too
  
  delta = range[:,1]-range[:,0]
  
  for i in xrange(len(data)):
    x = int(options.bins[0]*(data[i,0]-range[0,0])/delta[0])
    y = int(options.bins[1]*(data[i,1]-range[1,0])/delta[1])
    if x >=0 and x < options.bins[0] and y >= 0 and y < options.bins[1]: grid[x,y] += 1 if options.weight == None else data[i,2]
  
  if (range[2] == 0.0).all(): range[2] = [grid.min(),grid.max()]
  if (range[2] == 0.0).all():                                             # no data in grid?
    file['croak'].write('no data found on grid...\n')
    range[2,:] = numpy.array([0.0,1.0])                                   # making up arbitrary z range
  if options.type[2].lower() == 'log':
    grid = numpy.log(grid)
    range[2] = numpy.log(range[2])
    
  delta[2] = range[2,1]-range[2,0]

  
  i = 0
  for x in xrange(options.bins[0]):
    for y in xrange(options.bins[1]):
      result[i,:] = [range[0,0]+delta[0]/options.bins[0]*(x+0.5),
                     range[1,0]+delta[1]/options.bins[1]*(y+0.5),
                     min(1.0,max(0.0,(grid[x,y]-range[2,0])/delta[2]))]
      if options.type[0].lower() == 'log': result[i,0] = numpy.exp(result[i,0])
      if options.type[1].lower() == 'log': result[i,1] = numpy.exp(result[i,1])
      if options.invert:                   result[i,2] = 1.0-result[i,2]
      i += 1
        
  file['output'].write('1\thead\n')
  file['output'].write('bin_%s\tbin_%s\tz\n'%(headers[options.data[0]-1],headers[options.data[1]-1]))
  numpy.savetxt(file['output'],result)
  file['output'].close()                                                  # close output ASCII table
