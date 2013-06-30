#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math,numpy
import damask
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

synonyms = {
        'grid':   ['resolution'],
        'size':   ['dimension'],
          }
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
Smoothens out interface roughness by simulated curvature flow.
This is achieved by the diffusion of each initially sharply bounded grain volume within the periodic domain
up to a given distance 'd' voxels.
The final geometry is assembled by selecting at each voxel that grain index for which the concentration remains largest.
""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-d', '--distance', dest='d', type='float', \
                  help='diffusion distance in voxels [%default]', metavar='float')
parser.add_option('-b', '--black', dest='black', action='extend', type='string', \
                  help='indices of stationary microstructures', metavar='<LIST>')

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

  theTable = damask.ASCIItable(file['input'],file['output'],labels=False)
  theTable.head_read()

#--- interpret header ----------------------------------------------------------------------------
  info = {
          'grid':   numpy.zeros(3,'i'),
          'size':   numpy.zeros(3,'d'),
          'origin': numpy.zeros(3,'d'),
          'homogenization':  0,
          'microstructures': 0,
         }
  newInfo = {
          'microstructures': 0,
         }
  extra_header = []

  for header in theTable.info:
    headitems = map(str.lower,header.split())
    if len(headitems) == 0: continue
    for synonym,alternatives in synonyms.iteritems():
      if headitems[0] in alternatives: headitems[0] = synonym
    if headitems[0] in mappings.keys():
      if headitems[0] in identifiers.keys():
        for i in xrange(len(identifiers[headitems[0]])):
          info[headitems[0]][i] = \
            mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])
      else:
        info[headitems[0]] = mappings[headitems[0]](headitems[1])
    else:
      extra_header.append(header)

  file['croak'].write('grid     a b c:  %s\n'%(' x '.join(map(str,info['grid']))) + \
                      'size     x y z:  %s\n'%(' x '.join(map(str,info['size']))) + \
                      'origin   x y z:  %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization:  %i\n'%info['homogenization'] + \
                      'microstructures: %i\n'%info['microstructures'])

  if numpy.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    continue
  if numpy.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    continue

#--- read data ------------------------------------------------------------------------------------
  microstructure = numpy.zeros(info['grid'].prod(),'i')
  i = 0
  theTable.data_rewind()
  while theTable.data_read():
    items = theTable.data
    if len(items) > 2:
      if   items[1].lower() == 'of': items = [int(items[2])]*int(items[0])
      elif items[1].lower() == 'to': items = xrange(int(items[0]),1+int(items[2]))
      else:                          items = map(int,items)
    else:                            items = map(int,items)

    s = len(items)
    microstructure[i:i+s] = items
    i += s

#--- do work ------------------------------------------------------------------------------------
  microstructure = microstructure.reshape(info['grid'],order='F')

#--- initialize helper data -----------------------------------------------------------------------  
  winner      = numpy.zeros(info['grid'],'i')
  diffusedMax = numpy.zeros(info['grid'])

#--- diffuse each grain separately ----------------------------------------------------------------  
  for theGrain in xrange(1,1+numpy.amax(microstructure)):
    diffused = ndimage.filters.gaussian_filter((microstructure == theGrain).astype(float),\
                                               {True:0.0,False:options.d}[theGrain in options.black],\
                                               mode='wrap')
    winner = numpy.where(diffused > diffusedMax, theGrain, winner)
    diffusedMax = numpy.where(diffused > diffusedMax, diffused, diffusedMax)

  newInfo['microstructures'] = winner.max()

#--- report ---------------------------------------------------------------------------------------
  if (newInfo['microstructures'] != info['microstructures']):
    file['croak'].write('--> microstructures: %i\n'%newInfo['microstructures'])

#--- write header ---------------------------------------------------------------------------------
  theTable.labels_clear()
  theTable.info_clear()
  theTable.info_append(extra_header+[
    "$Id$",
    "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
    "size\tx %f\ty %f\tz %f"%(info['size'][0],info['size'][1],info['size'][2],),
    "origin\tx %f\ty %f\tz %f"%(info['origin'][0],info['origin'][1],info['origin'][2],),
    "homogenization\t%i"%info['homogenization'],
    "microstructures\t%i"%(newInfo['microstructures']),
    ])
  theTable.head_write()
  theTable.output_flush()
  
# --- write microstructure information ------------------------------------------------------------
  formatwidth = int(math.floor(math.log10(winner.max())+1))
  theTable.data = winner.reshape((info['grid'][0],info['grid'][1]*info['grid'][2]),order='F').transpose()
  theTable.data_writeArray('%%%ii'%(formatwidth))
    
#--- output finalization --------------------------------------------------------------------------
  if file['name'] != 'STDIN':
    file['input'].close()
    file['output'].close()
    os.rename(file['name']+'_tmp',file['name'])
