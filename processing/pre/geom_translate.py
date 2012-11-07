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
translate microstructure indices (shift or substitute) and/or geometry origin.
""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-o', '--origin', dest='origin', type='float', nargs = 3, \
                  help='offset from old to new origin of grid', metavar='<x y z>')
parser.add_option('-m', '--microstructure', dest='microstructure', type='int', \
                  help='offset from old to new microstructure indices', metavar='<int>')
parser.add_option('-s', '--substitute', action='extend', dest='substitute', type='string', \
                  help='substitutions of microstructure indices from,to,from,to,...', metavar='<list>')
parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                  help='output geom file with two-dimensional data arrangement')

parser.set_defaults(origin = [0.0,0.0,0.0])
parser.set_defaults(microstructure = 0)
parser.set_defaults(substitute = [])
parser.set_defaults(twoD = False)

(options, filenames) = parser.parse_args()

sub = {}
for i in xrange(len(options.substitute)/2):                                               # split substitution list into "from" -> "to"
  sub[int(options.substitute[i*2])] = int(options.substitute[i*2+1])

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

  info = {'resolution': numpy.array([0,0,0]),
          'dimension':  numpy.array([0.0,0.0,0.0]),
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

  file['croak'].write('resolution:     %s\n'%(' x '.join(map(str,info['resolution']))) + \
                      'dimension:      %s\n'%(' x '.join(map(str,info['dimension']))) + \
                      'origin:         %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization: %i\n'%info['homogenization'])

  new_header.append("resolution\ta %i\tb %i\tc %i\n"%( 
    info['resolution'][0],
    info['resolution'][1],
    info['resolution'][2],))
  new_header.append("dimension\tx %f\ty %f\tz %f\n"%(
    info['dimension'][0],
    info['dimension'][1],
    info['dimension'][2]))
  new_header.append("origin\tx %f\ty %f\tz %f\n"%(
    info['origin'][0]+options.origin[0],
    info['origin'][1]+options.origin[1],
    info['origin'][2]+options.origin[2]))
  new_header.append("homogenization\t%i\n"%info['homogenization'])

# ------------------------------------------ assemble header ---------------------------------------  

  output  = '%i\theader\n'%(len(new_header))
  output += ''.join(new_header)
  file['output'].write(output)
  
# ------------------------------------------ process input ---------------------------------------  

  N = info['resolution'][0]*info['resolution'][1]*info['resolution'][2]
  microstructure = numpy.zeros(N,'i')

  i = 0
  for line in content:  
    d = map(int,line.split())
    s = len(d)
    microstructure[i:i+s] = d                                                       # read microstructure indices
    i += s

  for i in xrange(N):
    if microstructure[i] in sub: microstructure[i] = sub[microstructure[i]]         # substitute microstructure indices
    
  microstructure += options.microstructure                                          # shift microstructure indices

  formatwidth = int(math.floor(math.log10(microstructure.max())+1))
  i = 0
  for z in xrange(info['resolution'][2]):
    for y in xrange(info['resolution'][1]):
      output = {True:' ',False:'\n'}[options.twoD].join(map(lambda x: ('%%%ii'%formatwidth)%x, microstructure[i:i+info['resolution'][0]])) + '\n'
      file['output'].write(output)
      i += info['resolution'][0]
  

# ------------------------------------------ output finalization ---------------------------------------  

  if file['name'] != 'STDIN':
    file['output'].close()
    os.rename(file['name']+'_tmp',file['name'])
