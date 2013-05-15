#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math,numpy
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP

#--------------------------------------------------------------------------------------------------
class extendedOption(Option):
#--------------------------------------------------------------------------------------------------
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

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
identifiers = {
        'grid':    ['a','b','c'],
        'size':    ['x','y','z'],
        'origin':  ['x','y','z'],
          }
mappings = {
        'grid':           lambda x: int(x),
        'size':           lambda x: float(x),
        'origin':         lambda x: float(x),
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
for i in xrange(len(options.substitute)/2):                                                         # split substitution list into "from" -> "to"
  sub[int(options.substitute[i*2])] = int(options.substitute[i*2+1])

#--- setup file handles ---------------------------------------------------------------------------  
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

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write(file['name']+'\n')

  firstline = file['input'].readline()
  m = re.search('(\d+)\s*head', firstline.lower())
  if m:
    headerlines = int(m.group(1))
    headers  = [file['input'].readline() for i in range(headerlines)]
  else:
    headerlines = 1
    headers = firstline

  content = file['input'].readlines()
  file['input'].close()

#--- interprete header ----------------------------------------------------------------------------
  info = {
          'grid':    numpy.zeros(3,'i'),
          'size':    numpy.zeros(3,'d'),
          'origin':  numpy.zeros(3,'d'),
          'microstructures': 0,
          'homogenization':  0
         }
  newInfo = {
          'origin': numpy.zeros(3,'d'),
          'microstructures': 0,
         }
  new_header = []
  for header in headers:
    headitems = map(str.lower,header.split())
    if headitems[0] == 'resolution': headitems[0] = 'grid'
    if headitems[0] == 'dimension':  headitems[0] = 'size'
    if headitems[0] in mappings.keys():
      if headitems[0] in identifiers.keys():
        for i in xrange(len(identifiers[headitems[0]])):
          info[headitems[0]][i] = \
            mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])
      else:
        info[headitems[0]] = mappings[headitems[0]](headitems[1])
    else:
      new_header.append(header)


  file['croak'].write('grid     a b c:  %s\n'%(' x '.join(map(str,info['grid']))) + \
                      'size     x y z:  %s\n'%(' x '.join(map(str,info['size']))) + \
                      'origin   x y z:  %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization:  %i\n'%info['homogenization'] + \
                      'microstructures: %i\n\n'%info['microstructures'])

  if numpy.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    sys.exit()
  if numpy.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    sys.exit()

#--- process input --------------------------------------------------------------------------------  
  N = info['grid'].prod()
  microstructure = numpy.zeros(N,'i')

  i = 0
  for line in content:  
    d = map(int,line.split())
    s = len(d)
    microstructure[i:i+s] = d                                                                       # read microstructure indices
    i += s

  for i in xrange(N):
    if microstructure[i] in sub: microstructure[i] = sub[microstructure[i]]                         # substitute microstructure indices
    
  microstructure += options.microstructure                                                          # shift microstructure indices
  
  formatwidth = int(math.floor(math.log10(microstructure.max())+1))

#--- assemble header and report changes -----------------------------------------------------------
  newInfo['origin'] = info['origin'] + options.origin
  newInfo['microstructures'] = microstructure.max()

  if (any(newInfo['origin'] != info['origin'])):
    file['croak'].write('--> origin   x y z:  %s\n'%(' : '.join(map(str,newInfo['origin']))))
  if (newInfo['microstructures'] != info['microstructures']):
    file['croak'].write('--> microstructures: %i\n'%newInfo['microstructures'])

  new_header.append('$Id$\n')
  new_header.append("grid\ta %i\tb %i\tc %i\n"%(info['grid'][0],info['grid'][1],info['grid'][2],))
  new_header.append("size\tx %f\ty %f\tz %f\n"%(info['size'][0],info['size'][1],info['size'][2],))
  new_header.append("origin\tx %f\ty %f\tz %f\n"%(
                                 newInfo['origin'][0],newInfo['origin'][1],newInfo['origin'][2],))
  new_header.append("microstructures\t%i\n"%newInfo['microstructures'])
  new_header.append("homogenization\t%i\n"%info['homogenization'])
  file['output'].write('%i\theader\n'%(len(new_header))+''.join(new_header))

#--- write new data ------------------------------------------------------------------------------- 
  i = 0
  for z in xrange(info['grid'][2]):
    for y in xrange(info['grid'][1]):
      output = {True:' ',False:'\n'}[options.twoD].join(map(lambda x: ('%%%ii'%formatwidth)%x, 
                                                         microstructure[i:i+info['grid'][0]])) + '\n'
      file['output'].write(output)
      i += info['grid'][0]

#--- output finalization --------------------------------------------------------------------------  
  if file['name'] != 'STDIN':
    file['output'].close()
    os.rename(file['name']+'_tmp',file['name'])
