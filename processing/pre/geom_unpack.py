#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
synonyms = {
        'grid':   ['resolution'],
        'size':   ['dimension'],
          }
identifiers = {
        'grid':      ['a','b','c'],
        'size':      ['x','y','z'],
        'origin':    ['x','y','z'],
          }
mappings = {
        'grid':            lambda x: int(x),
        'size':            lambda x: float(x),
        'origin':          lambda x: float(x),
        'homogenization':  lambda x: int(x),
        'microstructures': lambda x: int(x),
          }

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Unpack geometry files containing ranges "a to b" and/or "n of x" multiples (exclusively in one line).
""", version = scriptID)

parser.add_option('-1', '--onedimensional', dest='oneD', action='store_true', \
                  help='output geom file with one-dimensional data arrangement [%default]')

parser.set_defaults(oneD = False)

(options, filenames) = parser.parse_args()

# ------------------------------------------ setup file handles ------------------------------------
files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

# ------------------------------------------ loop over input files ---------------------------------
for file in files:
  file['croak'].write('\033[1m' + scriptName + '\033[0m: ' + (file['name'] if file['name'] != 'STDIN' else '') + '\n')

  table = damask.ASCIItable(file['input'],file['output'],labels = False,buffered = False)           # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info


#--- interpret header ----------------------------------------------------------------------------
  info = {
          'grid':    np.zeros(3,'i'),
          'size':    np.zeros(3,'d'),
          'origin':  np.zeros(3,'d'),
          'homogenization':  0,
          'microstructures': 0,
         }
  extra_header = []

  for header in table.info:
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

  if np.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    continue
  if np.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    continue

#--- read data ------------------------------------------------------------------------------------
  microstructure = np.zeros(info['grid'].prod(),'i')
  i = 0

  while table.data_read():                                                                         # read next data line of ASCII table
    items = table.data
    if len(items) > 2:
      if   items[1].lower() == 'of': items = [int(items[2])]*int(items[0])
      elif items[1].lower() == 'to': items = xrange(int(items[0]),1+int(items[2]))
      else:                          items = map(int,items)
    else:                            items = map(int,items)

    s = len(items)
    microstructure[i:i+s] = items
    i += s

#--- write header ---------------------------------------------------------------------------------
  table.labels_clear()
  table.info_clear()
  table.info_append(extra_header+[
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
    "size\tx %e\ty %e\tz %e"%(info['size'][0],info['size'][1],info['size'][2],),
    "origin\tx %e\ty %e\tz %e"%(info['origin'][0],info['origin'][1],info['origin'][2],),
    "homogenization\t%i"%info['homogenization'],
    "microstructures\t%i"%(info['microstructures']),
    ])
  table.head_write()
  
# --- write microstructure information ------------------------------------------------------------
  formatwidth = int(math.floor(math.log10(microstructure.max())+1))
  if options.oneD:
    table.data = microstructure
  else:
    table.data = microstructure.reshape((info['grid'][0],info['grid'][1]*info['grid'][2]),order='F').transpose()
  table.data_writeArray('%%%ii'%(formatwidth),delimiter = ' ')
    
#--- output finalization --------------------------------------------------------------------------
  if file['name'] != 'STDIN':
    table.input_close()                                                                            # close input ASCII table
    table.output_close()                                                                           # close input ASCII table
    os.rename(file['name']+'_tmp',file['name'])                                                    # overwrite old one with tmp new
