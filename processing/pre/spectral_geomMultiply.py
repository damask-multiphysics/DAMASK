#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math,numpy
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP


# -----------------------------
class extendedOption(Option):
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


# ----------------------- MAIN -------------------------------

identifiers = {
        'resolution': ['a','b','c'],
        'dimension':  ['x','y','z'],
          }
mappings = {
        'resolution': lambda x: int(x),
        'dimension':  lambda x: float(x),
          }


parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
Multiplies a spectral geometry description independently in x, y, and z direction at fixed dimension.
""" + string.replace('$Id: spectral_geomCrop.py 1449 2012-04-24 11:31:18Z MPIE\t.elachkar $','\n','\\n')
)

parser.add_option('-f', '--factor', dest='factor', type='int', nargs = 3, \
                  help='multiplication factors of a,b,c resolution')
parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                  help='output geom file with two-dimensional data arrangement')
parser.add_option('-k', '--keep', dest='keepDimension', action='store_true', \
                  help='keep dimension as before')

parser.set_defaults(factor = [1,1,1])
parser.set_defaults(twoD = False)
parser.set_defaults(keepDimension = False)

(options, filenames) = parser.parse_args()

prefix = 'mult%ix%ix%i'%(options.factor[0],options.factor[1],options.factor[2])

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(prefix+'_'+name,'w')})


# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': print file['name']

  #  get labels by either read the first row, or - if keyword header is present - the last line of the header

  firstline = file['input'].readline()
  m = re.search('(\d+)\s*head', firstline.lower())
  if m:
    headerlines = int(m.group(1))
    headers  = [firstline]+[file['input'].readline() for i in range(headerlines)]
  else:
    headerlines = 1
    headers = firstline

  content = file['input'].readlines()
  file['input'].close()

  resolution = [0,0,0]
  dimension  = [0.0,0.0,0.0]
  resolution_blown = [0,0,0]
  dimension_blown  = [0.0,0.0,0.0]
  new_header = []
  for header in headers:
    headitems = header.split()
    if headitems[0] == 'resolution':         # located resolution entry
      for i in xrange(3):
        resolution[i] = mappings['resolution'](headitems[headitems.index(identifiers['resolution'][i])+1])
        resolution_blown[i] = resolution[i]*options.factor[i]
      header = "resolution\ta %i\tb %i\tc %i\n"%(resolution_blown[0],resolution_blown[1],resolution_blown[2])
    if headitems[0] == 'dimension':          # located dimension entry
      for i in xrange(3):
        dimension[i] = mappings['dimension'](headitems[headitems.index(identifiers['dimension'][i])+1])
        if options.keepDimension: dimension_blown[i] = dimension[i]
        else:                     dimension_blown[i] = dimension[i]*options.factor[i]
      header = "dimension\tx %f\ty %f\tz %f\n"%(dimension_blown[0],dimension_blown[1],dimension_blown[2])
    
    new_header.append(header)
    
  if resolution == [0,0,0]:
    print 'no resolution info found.'
    sys.exit(1)
  if dimension == [0.0,0.0,0.0]:
    print 'no dimension info found.'
    sys.exit(1)

  if file['name'] != 'STDIN':
    print 'resolution: %s'%(' x '.join(map(str,resolution_blown)))
    print 'dimension:  %s'%(' x '.join(map(str,dimension_blown)))
    
  microstructure = numpy.zeros(resolution,'i')
  i = 0
  for line in content:  
    for item in map(int,line.split()):
      microstructure[i%resolution[0],(i/resolution[0])%resolution[1],i/resolution[0]/resolution[1]] = item
      i += 1
  
  formatwidth = int(math.floor(math.log10(microstructure.max())))
  

            
# ------------------------------------------ assemble header ---------------------------------------  

  output = ''.join(new_header)

# ------------------------------------- regenerate texture information ----------------------------------  

  for z in xrange(resolution[2]):
    for c in xrange(options.factor[2]):
      for y in xrange(resolution[1]):
        for b in xrange(options.factor[1]):
          for x in xrange(resolution[0]):
            for a in xrange(options.factor[0]):
              output += ('%%%ii'%formatwidth)%microstructure[x,y,z] + {True:' ',False:'\n'}[options.twoD]
          output += {True:'\n',False:''}[options.twoD]
    
# ------------------------------------------ output result ---------------------------------------  

  file['output'].write(output)

  file['input'].close()
  if file['name'] != 'STDIN':
    file['output'].close()
    