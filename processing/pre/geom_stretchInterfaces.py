#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math,numpy
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP
from scipy import ndimage

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
        'grid':           lambda x: int(x),
        'size':           lambda x: float(x),
        'origin':         lambda x: float(x),
        'homogenization': lambda x: int(x),
          }

parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
Smoothens out interface roughness by simulated curvature flow.
This is achieved by the diffusion of each initially sharply bounded grain volume within the periodic domain for a given time,
i.e. up to a diffusion distance of sqrt(t) voxels.
The final geometry is assembled by selecting at each voxel that grain index for which the concentration is largest.
""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-t', '--time', dest='t', type='int', \
                  help='time for curvature flow [%default]')
parser.add_option('-b', '--black', dest='black', action='extend', type='string', \
                  help='indices of stationary microstructures', metavar='<LIST>')
parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                  help='output geom file with two-dimensional data arrangement [%default]')

parser.set_defaults(t = 1)
parser.set_defaults(black = [])
parser.set_defaults(twoD = False)

(options, filenames) = parser.parse_args()

options.black = map(int,options.black)

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
          'homogenization':  0,
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
    items = line.split()
    if len(items) > 2:
      if   items[1].lower() == 'of': items = [int(items[2])]*int(items[0])
      elif items[1].lower() == 'to': items = xrange(int(items[0]),1+int(items[2]))
      else:                            items = map(int,items)
    else:                              items = map(int,items)

    for item in items:
      microstructure[i%info['grid'][0],
                    (i/info['grid'][0])%info['grid'][1],
                     i/info['grid'][0] /info['grid'][1]] = item
      i += 1

#--- initialize helper data -----------------------------------------------------------------------  
  winner      = numpy.zeros(info['grid'],'i')
  diffusedMax = numpy.zeros(info['grid'])

#--- diffuse each grain separately ----------------------------------------------------------------  
  for theGrain in xrange(1,1+numpy.amax(microstructure)):
    diffused = ndimage.filters.gaussian_filter((microstructure == theGrain).astype(float),\
                                               {True:0.0,False:numpy.sqrt(options.t)}[theGrain in options.black],\
                                               mode='wrap')
    winner = numpy.where(diffused > diffusedMax, theGrain, winner)
    diffusedMax = numpy.where(diffused > diffusedMax, diffused, diffusedMax)

  microstructure = winner


# --- assemble header -----------------------------------------------------------------------------
  formatwidth = int(math.floor(math.log10(microstructure.max())+1))

  new_header.append('$Id$\n')
  new_header.append("grid\ta %i\tb %i\tc %i\n"%(
                     info['grid'][0],info['grid'][1],info['grid'][2]))
  new_header.append("size\tx %f\ty %f\tz %f\n"%(
                     info['size'][0],info['size'][1],info['size'][2]))
  new_header.append("origin\tx %f\ty %f\tz %f\n"%(
                     info['origin'][0],info['origin'][1],info['origin'][2]))
  new_header.append("homogenization\t%i\n"%info['homogenization'])
  new_header.append("microstructures\t%i\n"%info['microstructures'])
  file['output'].write('%i\theader\n'%(len(new_header))+''.join(new_header))

# --- write microstructure information ------------------------------------------------------------
  for z in xrange(info['grid'][2]):
    for y in xrange(info['grid'][1]):
      file['output'].write({True:' ',False:'\n'}[options.twoD].join(map(lambda x: \
                                    ('%%%ii'%formatwidth)%x, microstructure[:,y,z])) + '\n')

#--- output finalization --------------------------------------------------------------------------
  if file['name'] != 'STDIN':
    file['output'].close()
    os.rename(file['name']+'_tmp',file['name'])
