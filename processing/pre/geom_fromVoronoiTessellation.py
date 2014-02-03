#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,sys,math,string,re,numpy
import damask 
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
  

def meshgrid2(*arrs):
  '''
  code inspired by http://stackoverflow.com/questions/1827489/numpy-meshgrid-in-3d
  '''
  arrs = tuple(reversed(arrs))
  arrs = tuple(arrs)
  lens = numpy.array(map(len, arrs))
  dim = len(arrs)
  ans = []
  for i, arr in enumerate(arrs):
     slc = numpy.ones(dim,'i')
     slc[i] = lens[i]
     arr2 = numpy.asarray(arr).reshape(slc)
     for j, sz in enumerate(lens):
         if j != i:
             arr2 = arr2.repeat(sz, axis=j)
   
     ans.insert(0,arr2)
  return tuple(ans)


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
synonyms = {
        'grid':              ['resolution'],
        'size':              ['dimension'],
        'microstructures':   ['grains'],
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
Generate geometry description and material configuration by standard Voronoi tessellation of given seeds file.
""" + string.replace(scriptID,'\n','\\n')
)

parser.add_option('-g', '--grid', dest='grid', type='int', nargs = 3, metavar = 'int int int', \
                  help='a,b,c grid of hexahedral box [from seeds file]')
parser.add_option('-s', '--size', dest='size', type='float', nargs = 3, metavar = 'float float float', \
                  help='x,y,z size of hexahedral box [1.0 along largest grid point number]')
parser.add_option('--homogenization', dest='homogenization', type='int', metavar = 'int', \
                  help='homogenization index to be used [%default]')
parser.add_option('--phase', dest='phase', type='int', metavar = 'int', \
                  help='phase index to be used [%default]')
parser.add_option('--crystallite', dest='crystallite', type='int', metavar = 'int', \
                  help='crystallite index to be used [%default]')
parser.add_option('-c', '--configuration', dest='config', action='store_true', \
                  help='output material configuration [%default]')
                                   
parser.set_defaults(grid = [0,0,0])
parser.set_defaults(size = [0.0,0.0,0.0])
parser.set_defaults(homogenization = 1)
parser.set_defaults(phase          = 1)
parser.set_defaults(crystallite    = 1)
parser.set_defaults(config = False)

(options,filenames) = parser.parse_args()

#--- setup file handles ---------------------------------------------------------------------------  
files = []
if filenames == []:
  files.append({'name':   'STDIN',
                'input':  sys.stdin,
                'output': sys.stdout,
                'croak':  sys.stderr,
               })
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':   name,
                    'input':  open(name),
                    'output': open(name+'_tmp','w'),
                    'croak':  sys.stdout,
                    })


#--- loop over input files ------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  theTable = damask.ASCIItable(file['input'],file['output'],buffered = False)
  theTable.head_read()

  labels = ['x','y','z']
  index = 0

  hasEulers = numpy.all(theTable.labels_index(['phi1','Phi','phi2'])) != -1
  hasGrains = theTable.labels_index('microstructure') != -1

  if hasEulers:
    labels += ['phi1','Phi','phi2']
    index += 3

  eulerCol = index

  if hasGrains:
    labels += ['microstructure']
    index += 1

  grainCol = index

  theTable.data_readArray(labels)
  coords = theTable.data[:,0:3]
  eulers = theTable.data[:,eulerCol:eulerCol+3] if hasEulers else numpy.zeros(3*len(coords))
  grain = theTable.data[:,grainCol] if hasGrains else 1+numpy.arange(len(eulers))
  grainIDs = numpy.unique(grain).astype('i')


#--- interpret header ----------------------------------------------------------------------------
  info = {
          'grid':    numpy.zeros(3,'i'),
          'size':    numpy.array(options.size),
          'origin':  numpy.zeros(3,'d'),
          'microstructures':  0,
          'homogenization': options.homogenization,
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

  if info['microstructures'] != len(grainIDs):
    file['croak'].write('grain data not matching grain count (%i)...\n'%(len(grainIDs)))
    info['microstructures'] = len(grainIDs)
  
  if 0 not in options.grid:                                                                         # user-specified grid
    info['grid'] = numpy.array(options.grid)

  for i in xrange(3):
    if info['size'][i] <= 0.0:                                                                      # any invalid size?
      info['size'][i] = float(info['grid'][i])/max(info['grid'])
      file['croak'].write('rescaling size %s...\n'%{0:'x',1:'y',2:'z'}[i])

  file['croak'].write('grains to map:  %i\n'%info['microstructures'] + \
                      'grid     a b c: %s\n'%(' x '.join(map(str,info['grid']))) + \
                      'size     x y z: %s\n'%(' x '.join(map(str,info['size']))) + \
                      'origin   x y z: %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization: %i\n'%info['homogenization'])
  
  if numpy.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    continue
  if numpy.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    continue
  if info['microstructures'] == 0:
    file['croak'].write('no grain info found.\n')
    continue

#--- prepare data ---------------------------------------------------------------------------------
  coords = (coords*info['size']).transpose()
  eulers = eulers.transpose()

#--- switch according to task ---------------------------------------------------------------------
  if options.config:                                                                                # write config file
    formatwidth = 1+int(math.log10(info['microstructures']))
    file['output'].write('<microstructure>\n')
    for i in grainIDs:
      file['output'].write('\n[Grain%s]\n'%(str(i).zfill(formatwidth)) + \
                           'crystallite %i\n'%options.crystallite + \
                           '(constituent)\tphase %i\ttexture %s\tfraction 1.0\n'%(options.phase,str(i).rjust(formatwidth)))
  
    file['output'].write('\n<texture>\n')
    for i in grainIDs:
      eulerID = numpy.nonzero(grain == i)[0][0]                                                     # find first occurrence of this grain id
      file['output'].write('\n[Grain%s]\n'%(str(i).zfill(formatwidth)) + \
                           '(gauss)\tphi1 %g\tPhi %g\tphi2 %g\tscatter 0.0\tfraction 1.0\n'%(eulers[0,eulerID],
                                                                                             eulers[1,eulerID],
                                                                                             eulers[2,eulerID]))

  else:                                                                                             # write geometry file
    x = (numpy.arange(info['grid'][0])+0.5)*info['size'][0]/info['grid'][0]
    y = (numpy.arange(info['grid'][1])+0.5)*info['size'][1]/info['grid'][1]
    z = (numpy.arange(info['grid'][2])+0.5)*info['size'][2]/info['grid'][2]
    undeformed = numpy.vstack(map(numpy.ravel, meshgrid2(x, y, z)))

    file['croak'].write('tessellating...\n')
    indices = damask.core.math.periodicNearestNeighbor(\
              info['size'],\
              numpy.eye(3),\
              undeformed,coords)//3**3 + 1                                                          # floor division to kill periodic images
    indices = grain[indices-1]

    newInfo['microstructures'] = info['microstructures']
    for i in grainIDs:
      if i not in indices: newInfo['microstructures'] -= 1
    file['croak'].write({True:'all',False:'only'}[newInfo['microstructures']  == info['microstructures'] ] +
                        ' %i'%newInfo['microstructures'] + 
                        {True:'',False:' out of %i'%info['microstructures']}[newInfo['microstructures'] == info['microstructures']] +
                        ' grains mapped.\n')

#--- write header ---------------------------------------------------------------------------------
    theTable.labels_clear()
    theTable.info_clear()
    theTable.info_append(extra_header+[
      scriptID + ' ' + ' '.join(sys.argv[1:]),
      "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
      "size\tx %f\ty %f\tz %f"%(info['size'][0],info['size'][1],info['size'][2],),
      "origin\tx %f\ty %f\tz %f"%(info['origin'][0],info['origin'][1],info['origin'][2],),
      "homogenization\t%i"%info['homogenization'],
      "microstructures\t%i"%(newInfo['microstructures']),
      ])
    theTable.head_write()
    
# --- write microstructure information ------------------------------------------------------------
    formatwidth = 1+int(math.log10(newInfo['microstructures']))
    theTable.data = indices.reshape(info['grid'][1]*info['grid'][2],info['grid'][0])
    theTable.data_writeArray('%%%ii'%(formatwidth),delimiter=' ')
    
#--- output finalization --------------------------------------------------------------------------

  if file['name'] != 'STDIN':
    file['input'].close()
    file['output'].close()
    os.rename(file['name']+'_tmp',os.path.splitext(file['name'])[0] + \
                                  {True: '_material.config',
                                   False:'.geom'}[options.config])
