#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math,numpy
import damask
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP

scriptID = '$Id: geom_addPrimitive.py 3412 2014-08-22 16:58:53Z MPIE\m.diehl $'
scriptName = scriptID.split()[1]

oversampling = 2.

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
Positions a geometric object within the (three-dimensional) canvas of a spectral geometry description.
Depending on the sign of the dimension parameters, these objects can be boxes, cylinders, or ellipsoids.
""" + string.replace(scriptID,'\n','\\n')
)

parser.add_option('-o', '--origin', 
                  '-c', '--center',     dest='center', type='int', nargs = 3, metavar=' '.join(['int']*3),
                  help='a,b,c origin of primitive %default')
parser.add_option('-d', '--dimension',  dest='dimension', type='int', nargs = 3, metavar=' '.join(['int']*3),
                  help='a,b,c extension of hexahedral box; negative values are diameters')
parser.add_option('-f', '--fill',       dest='fill', type='int', metavar = 'int',
                  help='grain index to fill primitive. "0" selects maximum microstructure index + 1 [%default]')
parser.add_option('-q', '--quaternion', dest='quaternion', type='float', nargs = 4, metavar=' '.join(['float']*4),
                  help = 'rotation of primitive as quaternion')
parser.add_option('-a', '--angleaxis',  dest='angleaxis', nargs = 4, metavar=' '.join(['float']*4),
                  help = 'rotation of primitive as angle and axis')
parser.add_option(     '--degrees',     dest='degrees', action='store_true',
                  help = 'angle is given in degrees [%default]')

parser.set_defaults(center = [0,0,0],
                    fill = 0,
                    quaternion = [],
                    angleaxis = [],
                    degrees = False,
                   )

(options, filenames) = parser.parse_args()

if options.angleaxis != []:
  options.angleaxis = map(float,options.angleaxis)
  rotation = damask.Quaternion().fromAngleAxis(numpy.radians(options.angleaxis[0]) if options.degrees else options.angleaxis[0],
                                               options.angleaxis[1:4]).conjugated()
elif options.quaternion != []:
  options.rotation = map(float,options.rotation)
  rotation = damask.Quaternion(options.quaternion).conjugated()
else:
  rotation = damask.Quaternion().conjugated()

options.center = numpy.array(options.center)
invRotation = rotation.conjugated()                                                             # rotation of gridpos into primitive coordinate system

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

  table = damask.ASCIItable(file['input'],file['output'],labels = False)
  table.head_read()

#--- interpret header ----------------------------------------------------------------------------
  info = {
          'grid':    numpy.zeros(3,'i'),
          'size':    numpy.zeros(3,'d'),
          'origin':  numpy.zeros(3,'d'),
          'homogenization':  0,
          'microstructures': 0,
         }
  newInfo = {
          'grid':    numpy.zeros(3,'i'),
          'origin':  numpy.zeros(3,'d'),
          'microstructures': 0,
         }
  extra_header = []

  for header in table.info:
    headitems = map(str.lower,header.split())
    if len(headitems) == 0: continue                                                              # skip blank lines
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
  microstructure = numpy.zeros(info['grid'].prod(),'i')                                            # initialize as flat array
  i = 0

  while table.data_read():
    items = table.data
    if len(items) > 2:
      if   items[1].lower() == 'of': items = [int(items[2])]*int(items[0])
      elif items[1].lower() == 'to': items = xrange(int(items[0]),1+int(items[2]))
      else:                            items = map(int,items)
    else:                              items = map(int,items)

    s = len(items)
    microstructure[i:i+s] = items
    i += s

#--- do work ------------------------------------------------------------------------------------

  if options.fill == 0:
    options.fill = microstructure.max()+1

  microstructure = microstructure.reshape(info['grid'],order='F')

  if options.dimension != None:
    mask = (numpy.array(options.dimension) < 0).astype(float)                                       # zero where positive dimension, otherwise one
    dim = abs(numpy.array(options.dimension))                                                       # dimensions of primitive body
    pos = numpy.zeros(3,dtype='float')
#    hiresPrimitive = numpy.zeros((2*dim[0],2*dim[1],2*dim[2],3))                                   # primitive discretized at twice the grid resolution
    for     i,pos[0] in enumerate(numpy.arange(-dim[0]/oversampling,(dim[0]+1)/oversampling,1./oversampling)):
      for   j,pos[1] in enumerate(numpy.arange(-dim[1]/oversampling,(dim[1]+1)/oversampling,1./oversampling)):
        for k,pos[2] in enumerate(numpy.arange(-dim[2]/oversampling,(dim[2]+1)/oversampling,1./oversampling)):
          gridpos = numpy.floor(rotation*pos)                                                       # rotate and lock into spacial grid
          primPos = invRotation*gridpos                                                             # rotate back to primitive coordinate system
          if numpy.dot(mask*primPos/dim,mask*primPos/dim) <= 0.25 and \
             numpy.all(abs((1.-mask)*primPos/dim) <= 0.5):                                          # inside ellipsoid and inside box
             microstructure[(gridpos[0]+options.center[0])%info['grid'][0],
                            (gridpos[1]+options.center[1])%info['grid'][1],
                            (gridpos[2]+options.center[2])%info['grid'][2]] = options.fill          # assign microstructure index

  newInfo['microstructures'] = microstructure.max()


#--- report ---------------------------------------------------------------------------------------
  if (newInfo['microstructures'] != info['microstructures']):
    file['croak'].write('--> microstructures: %i\n'%newInfo['microstructures'])

#--- write header ---------------------------------------------------------------------------------
  table.labels_clear()
  table.info_clear()
  table.info_append(extra_header+[
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
    "size\tx %f\ty %f\tz %f"%(info['size'][0],info['size'][1],info['size'][2],),
    "origin\tx %f\ty %f\tz %f"%(info['origin'][0],info['origin'][1],info['origin'][2],),
    "homogenization\t%i"%info['homogenization'],
    "microstructures\t%i"%(newInfo['microstructures']),
    ])
  table.head_write()
  table.output_flush()
    
# --- write microstructure information ------------------------------------------------------------
  formatwidth = int(math.floor(math.log10(microstructure.max())+1))
  table.data = microstructure.reshape((info['grid'][0],info['grid'][1]*info['grid'][2]),order='F').transpose()
  table.data_writeArray('%%%ii'%(formatwidth),delimiter=' ')
    
#--- output finalization --------------------------------------------------------------------------
  if file['name'] != 'STDIN':
    table.input_close()  
    table.output_close()  
    os.rename(file['name']+'_tmp',file['name'])
