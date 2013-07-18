#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math,numpy
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP

scriptID = '$Id$'
scriptName = scriptID.split()[1]

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
        'grid':   ['a','b','c'],
        'size':   ['x','y','z'],
        'origin': ['x','y','z'],
          }
mappings = {
        'grid':            lambda x: int(x),
        'size':            lambda x: float(x),
        'origin':          lambda x: float(x),
        'homogenization':  lambda x: int(x),
        'microstructures': lambda x: int(x),
          }

parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
Offset microstructure index for points which see a microstructure different from themselves within a given (cubic) vicinity,
i.e. within the region close to a grain/phase boundary.
""" + string.replace(scriptID,'\n','\\n')
)

parser.add_option('-v', '--vicinity', dest='vicinity', type='int', \
                  help='voxel distance checked for presence of other microstructure [%default]')
parser.add_option('-m', '--microstructureoffset', dest='offset', type='int', \
                  help='offset (positive or negative) for tagged microstructure. '+
                       '"0" selects maximum microstructure index [%default]')
parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                  help='output geom file with two-dimensional data arrangement')

parser.set_defaults(vicinity = 1)
parser.set_defaults(offset   = 0)
parser.set_defaults(twoD = False)

(options, filenames) = parser.parse_args()

#--- setup file handles --------------------------------------------------------------------------   
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
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

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
        'grid':   numpy.zeros(3,'i'),
        'size':   numpy.zeros(3,'d'),
        'origin': numpy.zeros(3,'d'),
        'microstructures': 0,
        'homogenization':  0
       }
  newInfo = {
        'microstructures': 0,
       }

  new_header = []
  new_header.append(scriptID+'\n')
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
                      'microstructures: %i\n'%info['microstructures'])

  if numpy.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    sys.exit()
  if numpy.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    sys.exit()

#--- read data ------------------------------------------------------------------------------------
  microstructure = numpy.zeros(info['grid'],'i')
  i = 0
  for line in content:  
    for item in map(int,line.split()):
      microstructure[i%info['grid'][0],
                    (i/info['grid'][0])%info['grid'][1],
                     i/info['grid'][0] /info['grid'][1]] = item
      i += 1


  if options.offset == 0:
    options.offset = microstructure.max()
  formatwidth = 1+int(math.floor(math.log10(abs(microstructure.max()+options.offset))))

#--- search ---------------------------------------------------------------------------------------
  for x in xrange(info['grid'][0]):
    for y in xrange(info['grid'][1]):
      for z in xrange(info['grid'][2]):

        me = microstructure[x,y,z]
        breaker = False
        
        for dx in xrange(-options.vicinity,options.vicinity+1):
          for dy in xrange(-options.vicinity,options.vicinity+1):
            for dz in xrange(-options.vicinity,options.vicinity+1):

              they = microstructure[(x+dx)%info['grid'][0],(y+dy)%info['grid'][1],(z+dz)%info['grid'][2]]
              if they != me and they != me+options.offset:                                          # located alien microstructure in vicinity
                microstructure[x,y,z] += options.offset                                             # tag myself as close to aliens!
                breaker = True
                break

            if breaker: break

          if breaker: break

  newInfo['microstructures'] = microstructure.max()
  if (newInfo['microstructures'] != info['microstructures']):
    file['croak'].write('--> microstructures: %i\n'%newInfo['microstructures'])

# --- assemble header -----------------------------------------------------------------------------         
  new_header.append("grid\ta %i\tb %i\tc %i\n"%(info['grid'][0],info['grid'][1],info['grid'][2]))
  new_header.append("size\tx %f\ty %f\tz %f\n"%(info['size'][0],info['size'][1],info['size'][0]))
  new_header.append("origin\tx %f\ty %f\tz %f\n"%(info['origin'][0],info['origin'][1],info['origin'][2]))
  new_header.append("microstructures\t%i\n"%newInfo['microstructures'])
  new_header.append("homogenization\t%i\n"%info['homogenization'])
  file['output'].write('%i\theader\n'%(len(new_header))+''.join(new_header))

# --- write microstructure information ------------------------------------------------------------
  for z in xrange(info['grid'][2]):
    for y in xrange(info['grid'][1]):
      file['output'].write({True:' ',False:'\n'}[options.twoD].
                       join(map(lambda x: str(x).rjust(formatwidth), microstructure[:,y,z])) + '\n')
    
    file['output'].write('\n')

#--- output finalization --------------------------------------------------------------------------
  if file['name'] != 'STDIN':
    file['output'].close()
    os.rename(file['name']+'_tmp',file['name'])
