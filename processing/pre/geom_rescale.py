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
        'origin':     ['x','y','z'],
          }
mappings = {
        'resolution': lambda x: int(x),
        'dimension':  lambda x: float(x),
        'origin':     lambda x: float(x),
        'homogenization': lambda x: int(x),
          }


parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
Scales a geometry description independently in x, y, and z direction in terms of resolution and/or dimension.
""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-r', '--resolution', dest='resolution', type='int', nargs = 3, \
                  help='new resolution (a,b,c)')
parser.add_option('-d', '--dimension', dest='dimension', type='float', nargs = 3, \
                  help='new dimension (x,y,z)')
parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                  help='output geom file with two-dimensional data arrangement')

parser.set_defaults(resolution = [0,0,0])
parser.set_defaults(dimension  = [0.0,0.0,0.0])
parser.set_defaults(twoD = False)

(options, filenames) = parser.parse_args()

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN',
                'input':sys.stdin,
                'output':sys.stdout,
                'croak':sys.stderr,
               })
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name,
                    'input':open(name),
                    'output':open(name+'_tmp','w'),
                    'croak':sys.stdout,
                    })

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': file['croak'].write(file['name']+'\n')

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

  info = {'resolution': numpy.array(options.resolution),
          'dimension':  numpy.array(options.dimension),
          'origin':     numpy.array([0.0,0.0,0.0]),
          'homogenization': 1,
         }

  new_header = []
  for header in headers:
    headitems = map(str.lower,header.split())
    if headitems[0] in mappings.keys():
      if headitems[0] in identifiers.keys():
        for i in xrange(len(identifiers[headitems[0]])):
          info[headitems[0]][i] = \
            mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])
      else:
        info[headitems[0]] = mappings[headitems[0]](headitems[1])
    
  if numpy.all(info['resolution'] == 0):
    file['croak'].write('no resolution info found.\n')
    continue
  if numpy.all(info['dimension'] == 0.0):
    file['croak'].write('no dimension info found.\n')
    continue

  if options.resolution == [0,0,0]:
    options.resolution = info['resolution']
  if options.dimension == [0.0,0.0,0.0]:
    options.dimension = info['dimension']

  file['croak'].write('resolution:     %s\n'%(' x '.join(map(str,info['resolution']))) + \
                      'dimension:      %s\n'%(' x '.join(map(str,info['dimension']))) + \
                      'origin:         %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization: %i\n'%info['homogenization'])

  new_header.append("resolution\ta %i\tb %i\tc %i\n"%( 
    options.resolution[0],
    options.resolution[1],
    options.resolution[2],))
  new_header.append("dimension\tx %f\ty %f\tz %f\n"%(
    options.dimension[0],
    options.dimension[1],
    options.dimension[2],))
  new_header.append("origin\tx %f\ty %f\tz %f\n"%(
    info['origin'][0],
    info['origin'][1],
    info['origin'][2],))
  new_header.append("homogenization\t%i\n"%info['homogenization'])
    
  microstructure = numpy.zeros(info['resolution'],'i')
  i = 0
  for line in content:  
    for item in map(int,line.split()):
      microstructure[i%info['resolution'][0],
                    (i/info['resolution'][0])%info['resolution'][1],
                     i/info['resolution'][0] /info['resolution'][1]] = item
      i += 1
  
  formatwidth = 1+int(math.floor(math.log10(microstructure.max())))
           
# ------------------------------------------ assemble header ---------------------------------------  

  output  = '%i\theader\n'%(len(new_header))
  output += ''.join(new_header)

# ------------------------------------- regenerate texture information ----------------------------------  

  for c in xrange(options.resolution[2]):
    z = int(info['resolution'][2]*(c+0.5)/options.resolution[2])%info['resolution'][2]
    for b in xrange(options.resolution[1]):
      y = int(info['resolution'][1]*(b+0.5)/options.resolution[1])%info['resolution'][1]
      for a in xrange(options.resolution[0]):
        x = int(info['resolution'][0]*(a+0.5)/options.resolution[0])%info['resolution'][0]
        output += str(microstructure[x,y,z]).rjust(formatwidth) + {True:' ',False:'\n'}[options.twoD]
      output += {True:'\n',False:''}[options.twoD]
    
# ------------------------------------------ output result ---------------------------------------  

  file['output'].write(output)

  if file['name'] != 'STDIN':
    file['input'].close()
    file['output'].close()
    os.rename(file['name']+'_tmp',file['name'])
    