#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math,numpy,itertools
import damask
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP
from scipy import ndimage
from multiprocessing import Pool

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
""" + string.replace(scriptID,'\n','\\n')
)

parser.add_option('-d', '--distance', dest='d', type='int', metavar='int', \
                  help='diffusion distance in voxels [%default]')
parser.add_option('-N', '--smooth', dest='N', type='int', metavar='int', \
                 help='N for curvature flow [%default]')
parser.add_option('-r', '--renumber', dest='renumber', action='store_true', \
                  help='renumber microstructure indices from 1...N [%default]')

parser.set_defaults(d = 1)
parser.set_defaults(N = 1)
parser.set_defaults(renumber = False)
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

  theTable = damask.ASCIItable(file['input'],file['output'],labels = False,buffered = False)
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
  microstructure = numpy.zeros(numpy.prod([2 if i == 1 else i for i in info['grid']]),'i')              # 2D structures do not work
  i = 0

  while theTable.data_read():                                  # read next data line of ASCII table
    items = theTable.data
    if len(items) > 2:
      if   items[1].lower() == 'of': items = [int(items[2])]*int(items[0])
      elif items[1].lower() == 'to': items = xrange(int(items[0]),1+int(items[2]))
      else:                          items = map(int,items)
    else:                            items = map(int,items)

    s = len(items)
    microstructure[i:i+s] = items
    i += s

#--- reshape, if 2D make copy ---------------------------------------------------------------------
  nMicrostuctures = numpy.prod([2 if i == 1 else i for i in info['grid']])
  if nMicrostuctures > info['grid'].prod():
    microstructure[info['grid'].prod():nMicrostuctures] = microstructure[0:info['grid'].prod()]
  microstructure = microstructure.reshape([2 if i == 1 else i for i in info['grid']],order='F')

#--- domain decomposition -------------------------------------------------------------------------  
  stride = numpy.array([2 if i == 1 else i for i in info['grid']],'i')
  if numpy.any(numpy.floor(stride) != stride): 
    file['croak'].write('invalid domain decomposition.\n')
    continue

#--- initialize helper data -----------------------------------------------------------------------  
  window = 0
  X,Y,Z = numpy.mgrid[0:stride[0]+2*window,0:stride[1]+2*window,0:stride[2]+2*window]
  gauss = numpy.exp(-(X*X+Y*Y+Z*Z)/(2.0*options.d*options.d))/math.pow(2.0*numpy.pi*options.d*options.d,1.5)
  gauss[:,:,(stride[2]+2*window)/2::] = gauss[:,:,(stride[2]+2*window)/2-1::-1]
  gauss[:,(stride[1]+2*window)/2::,:] = gauss[:,(stride[1]+2*window)/2-1::-1,:]
  gauss[(stride[0]+2*window)/2::,:,:] = gauss[(stride[0]+2*window)/2-1::-1,:,:]
  gauss = numpy.fft.rfftn(gauss)
  for smoothIter in xrange(options.N):
    interfacialEnergy = lambda A,B: (A*B != 0)*(A != B)*1.0
    struc = ndimage.generate_binary_structure(3,1)
    microExt = numpy.zeros([microstructure.shape[0]+2,microstructure.shape[1]+2,microstructure.shape[2]+2])
    microExt[1:-1,1:-1,1:-1] = microstructure
    boundary = numpy.zeros(microstructure.shape)
    for i in range(3):
      for j in range(3):
        for k in range(3):
          boundary = numpy.maximum(boundary,
                                   interfacialEnergy(microstructure,microExt[i:microstructure.shape[0]+i,
                                                                             j:microstructure.shape[1]+j,
                                                                             k:microstructure.shape[2]+k]))
    index = ndimage.morphology.distance_transform_edt(boundary == 0.,return_distances = False,return_indices = True)
    boundary = numpy.fft.irfftn(numpy.fft.rfftn(numpy.where(ndimage.morphology.binary_dilation(boundary != 0.,
                                                                                               structure = struc,
                                                                                               iterations = 2*options.d-1),
                                                          boundary[index[0].flatten(),index[1].flatten(),index[2].flatten()].reshape(microstructure.shape),
                                                          0.))*gauss)                                          
    index = ndimage.morphology.distance_transform_edt(boundary >= 0.5,return_distances=False,return_indices=True)
    microstructure = microstructure[index[0].flatten(),index[1].flatten(),index[2].flatten()].reshape(microstructure.shape)

# --- renumber to sequence 1...Ngrains if requested ------------------------------------------------
#  http://stackoverflow.com/questions/10741346/numpy-frequency-counts-for-unique-values-in-an-array  
  if options.renumber:
    newID=0
    for microstructureID,count in enumerate(numpy.bincount(microstructure.reshape(info['grid'].prod()))):
      if count != 0:
        newID+=1
        microstructure=numpy.where(microstructure==microstructureID,newID,microstructure).reshape(microstructure.shape)
# --- assemble header -----------------------------------------------------------------------------
  newInfo['microstructures'] = microstructure[0:info['grid'][0],0:info['grid'][1],0:info['grid'][2]].max()

#--- report ---------------------------------------------------------------------------------------
  if (newInfo['microstructures'] != info['microstructures']):
    file['croak'].write('--> microstructures: %i\n'%newInfo['microstructures'])

#--- write header ---------------------------------------------------------------------------------
  theTable.labels_clear()
  theTable.info_clear()
  theTable.info_append(extra_header+[
    scriptID+ ' ' + ' '.join(sys.argv[1:]),
    "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
    "size\tx %f\ty %f\tz %f"%(info['size'][0],info['size'][1],info['size'][2],),
    "origin\tx %f\ty %f\tz %f"%(info['origin'][0],info['origin'][1],info['origin'][2],),
    "homogenization\t%i"%info['homogenization'],
    "microstructures\t%i"%(newInfo['microstructures']),
    ])
  theTable.head_write()
  
# --- write microstructure information ------------------------------------------------------------
  formatwidth = int(math.floor(math.log10(microstructure.max())+1))
  theTable.data = microstructure[0:info['grid'][0],0:info['grid'][1],0:info['grid'][2]].reshape(numpy.prod(info['grid']),order='F').transpose()
  theTable.data_writeArray('%%%ii'%(formatwidth),delimiter=' ')
    
#--- output finalization --------------------------------------------------------------------------
  if file['name'] != 'STDIN':
    file['input'].close()
    file['output'].close()
    os.rename(file['name']+'_tmp',file['name'])
