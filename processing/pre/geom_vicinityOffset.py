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
Offset microstructure index for points which see a microstructure different from themselves within a given (cubic) vicinity,
i.e. within the region close to a grain/phase boundary.
""" + string.replace('$Id: spectral_geomCheck 994 2011-09-05 13:38:10Z MPIE\p.eisenlohr $','\n','\\n')
)

parser.add_option('-v', '--vicinity', dest='vicinity', type='int', \
                  help='voxel distance checked for presence of other microstructure [%default]')
parser.add_option('-o', '--offset', dest='offset', type='int', \
                  help='integer offset for tagged microstructure [%default]')
parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                  help='output geom file with two-dimensional data arrangement')

parser.set_defaults(vicinity = 1)
parser.set_defaults(offset   = 0)
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
    
  microstructure = numpy.zeros(info['resolution'],'i')
  i = 0
  for line in content:  
    for item in map(int,line.split()):
      microstructure[i%info['resolution'][0],
                    (i/info['resolution'][0])%info['resolution'][1],
                     i/info['resolution'][0] /info['resolution'][1]] = item
      i += 1

  formatwidth = 1+int(math.floor(math.log10(abs(microstructure.max()+options.offset))))
  if options.offset == 0:
    options.offset = microstructure.max()
  file['croak'].write('offset:         %i\n'%options.offset)
  
  for x in xrange(info['resolution'][0]):
    for y in xrange(info['resolution'][1]):
      for z in xrange(info['resolution'][2]):

        me = microstructure[x,y,z]
        breaker = False
        
        for dx in xrange(-options.vicinity,options.vicinity+1):
          for dy in xrange(-options.vicinity,options.vicinity+1):
            for dz in xrange(-options.vicinity,options.vicinity+1):

              they = microstructure[(x+dx)%info['resolution'][0],(y+dy)%info['resolution'][1],(z+dz)%info['resolution'][2]]
              if they != me and they != me+options.offset:                    # located alien microstructure in vicinity
                microstructure[x,y,z] += options.offset                       # tag myself as close to aliens!
                breaker = True
                break

            if breaker: break

          if breaker: break

            
# ------------------------------------------ assemble header ---------------------------------------  

  output = ''.join(headers)

# ------------------------------------- regenerate texture information ----------------------------------  

  for z in xrange(info['resolution'][2]):
    for y in xrange(info['resolution'][1]):
      output += {True:' ',False:'\n'}[options.twoD].join(map(lambda x: str(x).rjust(formatwidth), microstructure[:,y,z])) + '\n'
    
    output += '\n'
    
# ------------------------------------------ output result ---------------------------------------  

  file['output'].write(output)

  if file['name'] != 'STDIN':
    file['output'].close
    os.rename(file['name']+'_tmp',file['name'])
