#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,vtk
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
identifiers = {
        'grid':    ['a','b','c'],
        'size':    ['x','y','z'],
        'origin':  ['x','y','z'],
          }
mappings = {
        'grid':            lambda x: int(x),
        'size':            lambda x: float(x),
        'origin':          lambda x: float(x),
        'homogenization':  lambda x: int(x),
        'microstructures': lambda x: int(x),
          }

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [geomfile[s]]', description = """
Produce ASCIItable of structure data from geom description

""", version = scriptID)

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
                    'croak':sys.stdout,
                    })

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN':
    file['input'] = open(file['name'])
    file['output'] = open(os.path.splitext(file['name'])[0]+'.txt','w')

  file['croak'].write('\033[1m' + scriptName + '\033[0m' + (': '+file['name'] if file['name'] != 'STDIN' else '') + '\n')

  theTable = damask.ASCIItable(file['input'],file['output'],labels = False)
  theTable.head_read()

#--- interpret header ----------------------------------------------------------------------------

  info = {
          'grid':   np.zeros(3,'i'),
          'size':   np.zeros(3,'d'),
          'origin': np.zeros(3,'d'),
          'homogenization':  0,
          'microstructures': 0,
         }

  for header in theTable.info:
    headitems = map(str.lower,header.split())
    if len(headitems) == 0: continue
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

  if np.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    continue
  if np.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    continue

# ------------------------------------------ assemble header ---------------------------------------

  theTable.labels_clear()
  theTable.labels_append(['%i_pos'%(i+1) for i in range(3)]+['microstructure'])

  theTable.head_write()

#--- generate grid --------------------------------------------------------------------------------

  xx = np.arange(float(info['grid'][0]))/info['grid'][0]*info['size'][0]+info['origin'][0]
  yy = np.arange(float(info['grid'][1]))/info['grid'][1]*info['size'][1]+info['origin'][1]
  zz = np.arange(float(info['grid'][2]))/info['grid'][2]*info['size'][2]+info['origin'][2]
  
#--- read microstructure information --------------------------------------------------------------

  i = 0
  outputAlive = True
  
  while outputAlive and theTable.data_read():
    items = theTable.data
    if len(items) > 2:
      if   items[1].lower() == 'of': items = [int(items[2])]*int(items[0])
      elif items[1].lower() == 'to': items = xrange(int(items[0]),1+int(items[2]))
      else:                          items = map(int,items)
    else:                            items = map(int,items)

    for item in items:
      theTable.data = [xx[ i%info['grid'][0]],
                       yy[(i/info['grid'][0])%info['grid'][1]],
                       zz[ i/info['grid'][0]/info['grid'][1]],
                       item]
      i += 1
      outputAlive = theTable.data_write()                                      # output processed line
      if not outputAlive: break

# ------------------------------------------ finalize output ---------------------------------------

  theTable.output_flush()                                                   # just in case of buffered ASCII table
  
  if file['name'] != 'STDIN':
    file['input'].close()                                                   # close input ASCII table
    file['output'].close()                                                  # close output ASCII table
