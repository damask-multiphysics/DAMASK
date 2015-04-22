#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

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


parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Create seed file taking microstructure indices from given geom file but excluding black-listed grains.

""", version = scriptID)

parser.add_option('-w','--white',   dest='whitelist', action='extend', \
                                    help='white list of grain IDs', metavar='<LIST>')
parser.add_option('-b','--black',   dest='blacklist', action='extend', \
                                    help='black list of grain IDs', metavar='<LIST>')

parser.set_defaults(whitelist = [])
parser.set_defaults(blacklist = [])

(options,filenames) = parser.parse_args()

options.whitelist = map(int,options.whitelist)
options.blacklist = map(int,options.blacklist)

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
                    'output':open(os.path.splitext(name)[0]+'.seeds','w'),
                    'croak':sys.stdout,
                    })

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  file['croak'].write('\033[1m' + scriptName + '\033[0m: ' + (file['name'] if file['name'] != 'STDIN' else '') + '\n')

  table = damask.ASCIItable(file['input'],file['output'],labels = False,buffered = False)
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
  if 'origin' not in info:
    info['origin'] = np.zeros(3)

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


# ------------------------------------------ assemble header ---------------------------------------  

  table.info = [
                   scriptID,
                   "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
                   "size\tx %i\ty %i\tz %i"%(info['size'][0],info['size'][1],info['size'][2],),
                   "origin\tx %i\ty %i\tz %i"%(info['origin'][0],info['origin'][1],info['origin'][2],),
                  ]
  table.labels_clear()
  table.labels_append(['1_coords','2_coords','3_coords','microstructure'])                    # implicitly switching label processing/writing on
  table.head_write()
  
#--- filtering of grain voxels ------------------------------------------------------------------------------------
  table.data_clear()
  i = 0
  outputDead = False
  coord = np.zeros(3,'d')
  for coord[2] in xrange(info['grid'][2]):
    for coord[1] in xrange(info['grid'][1]):
      for coord[0] in xrange(info['grid'][0]):
        if (options.whitelist == [] and options.blacklist == []) or \
           (options.whitelist != [] and microstructure[i]     in options.whitelist) or \
           (options.blacklist != [] and microstructure[i] not in options.blacklist):
          table.data = list((coord+0.5)/info['grid'])+[microstructure[i]]
          outputDead = not table.data_write()
        i += 1
        if outputDead: break
      if outputDead: break
    if outputDead: break

# ------------------------------------------ output result ---------------------------------------  

  outputDead or table.output_flush()                                        # just in case of buffered ASCII table

  table.input_close()                                                       # close input ASCII table
  if file['name'] != 'STDIN':
    table.output_close()                                                    # close output ASCII table
