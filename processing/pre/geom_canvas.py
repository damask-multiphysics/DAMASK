#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,string
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
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

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Changes the (three-dimensional) canvas of a spectral geometry description.

""", version = scriptID)

parser.add_option('-g', '--grid',   dest='grid', nargs = 3, metavar=' '.join(['string']*3),
                                    help='a,b,c grid of hexahedral box [unchanged]')
parser.add_option('-o', '--offset', dest='offset', type='int', nargs = 3, metavar=' '.join(['int']*3),
                                    help='a,b,c offset from old to new origin of grid %default')
parser.add_option('-f', '--fill',   dest='fill', type='int', metavar = 'int',
                                    help='(background) canvas grain index. "0" selects maximum microstructure index + 1 [%default]')

parser.set_defaults(grid = ['0','0','0'])
parser.set_defaults(offset = (0,0,0))
parser.set_defaults(fill = 0)

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
  file['croak'].write('\033[1m' + scriptName + '\033[0m: ' + (file['name'] if file['name'] != 'STDIN' else '') + '\n')

  table = damask.ASCIItable(file['input'],file['output'],labels = False)
  table.head_read()

#--- interpret header ----------------------------------------------------------------------------
  info = {
          'grid':    np.zeros(3,'i'),
          'size':    np.zeros(3,'d'),
          'origin':  np.zeros(3,'d'),
          'homogenization':  0,
          'microstructures': 0,
         }
  newInfo = {
          'grid':    np.zeros(3,'i'),
          'origin':  np.zeros(3,'d'),
          'microstructures': 0,
         }
  extra_header = []

  for header in table.info:
    headitems = map(str.lower,header.split())
    if len(headitems) == 0: continue                                                              # skip blank lines
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

  if np.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    continue
  if np.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    continue

#--- read data ------------------------------------------------------------------------------------
  microstructure = np.zeros(info['grid'].prod(),'i')                                            # initialize as flat array
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
  newInfo['grid'] = np.array(int(o*float(n.translate(None,'xX'))) if [n[-1].lower() == 'x' else int(n.translate(None,'xX'))}[n[-1].lower() == 'x'] for o,n in zip(info['grid'],options.grid)],'i')
  newInfo['grid'] = np.where(newInfo['grid'] <= 0  , info['grid'],newInfo['grid'])

  microstructure = microstructure.reshape(info['grid'],order='F')
  microstructure_cropped = np.zeros(newInfo['grid'],'i')
  microstructure_cropped.fill({True:options.fill,False:microstructure.max()+1}[options.fill>0])
  xindex = list(set(xrange(options.offset[0],options.offset[0]+newInfo['grid'][0])) & \
                                                               set(xrange(info['grid'][0])))
  yindex = list(set(xrange(options.offset[1],options.offset[1]+newInfo['grid'][1])) & \
                                                               set(xrange(info['grid'][1])))
  zindex = list(set(xrange(options.offset[2],options.offset[2]+newInfo['grid'][2])) & \
                                                               set(xrange(info['grid'][2])))
  translate_x = [i - options.offset[0] for i in xindex]
  translate_y = [i - options.offset[1] for i in yindex]
  translate_z = [i - options.offset[2] for i in zindex]
  microstructure_cropped[min(translate_x):(max(translate_x)+1),\
                         min(translate_y):(max(translate_y)+1),\
                         min(translate_z):(max(translate_z)+1)] \
        = microstructure[min(xindex):(max(xindex)+1),\
                         min(yindex):(max(yindex)+1),\
                         min(zindex):(max(zindex)+1)]

  newInfo['size']   = info['size']/info['grid']*newInfo['grid']
  newInfo['origin'] = info['origin']+info['size']/info['grid']*options.offset
  newInfo['microstructures'] = microstructure_cropped.max()

#--- report ---------------------------------------------------------------------------------------
  if (any(newInfo['grid'] != info['grid'])):
    file['croak'].write('--> grid     a b c:  %s\n'%(' x '.join(map(str,newInfo['grid']))))
  if (any(newInfo['size'] != info['size'])):
    file['croak'].write('--> size     x y z:  %s\n'%(' x '.join(map(str,newInfo['size']))))
  if (any(newInfo['origin'] != info['origin'])):
    file['croak'].write('--> origin   x y z:  %s\n'%(' : '.join(map(str,newInfo['origin']))))
  if (newInfo['microstructures'] != info['microstructures']):
    file['croak'].write('--> microstructures: %i\n'%newInfo['microstructures'])

  if np.any(newInfo['grid'] < 1):
    file['croak'].write('invalid new grid a b c.\n')
    continue
  if np.any(newInfo['size'] <= 0.0):
    file['croak'].write('invalid new size x y z.\n')
    continue

#--- write header ---------------------------------------------------------------------------------
  table.labels_clear()
  table.info_clear()
  table.info_append(extra_header+[
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta %i\tb %i\tc %i"%(newInfo['grid'][0],newInfo['grid'][1],newInfo['grid'][2],),
    "size\tx %f\ty %f\tz %f"%(newInfo['size'][0],newInfo['size'][1],newInfo['size'][2],),
    "origin\tx %f\ty %f\tz %f"%(newInfo['origin'][0],newInfo['origin'][1],newInfo['origin'][2],),
    "homogenization\t%i"%info['homogenization'],
    "microstructures\t%i"%(newInfo['microstructures']),
    ])
  table.head_write()
  table.output_flush()
    
# --- write microstructure information ------------------------------------------------------------
  formatwidth = int(math.floor(math.log10(microstructure_cropped.max())+1))
  table.data = microstructure_cropped.reshape((newInfo['grid'][0],newInfo['grid'][1]*newInfo['grid'][2]),order='F').transpose()
  table.data_writeArray('%%%ii'%(formatwidth),delimiter=' ')
    
#--- output finalization --------------------------------------------------------------------------
  if file['name'] != 'STDIN':
    table.input_close()  
    table.output_close()  
    os.rename(file['name']+'_tmp',file['name'])
