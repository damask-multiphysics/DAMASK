#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,random,math,numpy as np
import damask
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP

scriptID = '$Id$'
scriptName = scriptID.split()[1]

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
Create seeds file by poking at 45 degree through given geom file.
Mimics APS Beamline 34-ID-E DAXM poking.

""", version = scriptID)

parser.add_option('-N', '--points', dest='N', type='int', metavar='int', \
                  help='number of poking locations [%default]')
parser.add_option('-z', '--planes', dest='z', type='float', nargs = 2, metavar='float float', \
                  help='top and bottom z plane')
parser.add_option('-x', action='store_true', dest='x', \
                  help='poke 45 deg along x')
parser.add_option('-y', action='store_true', dest='y', \
                  help='poke 45 deg along y')

parser.set_defaults(x = False)
parser.set_defaults(y = False)
parser.set_defaults(N = 16)

(options,filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------
if filenames == []:
  filenames = ['STDIN']

for name in filenames:
  if name == 'STDIN':
    file = {'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m\n')
  else:
    if not os.path.exists(name): continue
    file = {'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  theTable = damask.ASCIItable(file['input'],file['output'],labels = False)
  theTable.head_read()


#--- interpret header ----------------------------------------------------------------------------
  info = {
          'grid':    np.zeros(3,'i'),
          'size':    np.zeros(3,'d'),
          'origin':  np.zeros(3,'d'),
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

  Nx = int(options.N/math.sqrt(options.N*info['size'][1]/info['size'][0]))
  Ny = int(options.N/math.sqrt(options.N*info['size'][0]/info['size'][1]))
  Nz = int((max(options.z)-min(options.z))/info['size'][2]*info['grid'][2])

  file['croak'].write('poking %i x %i x %i...\n'%(Nx,Ny,Nz))
  microstructure = microstructure.reshape(info['grid'],order='F')
  seeds = np.zeros((Nx*Ny*Nz,4),'d')
  grid = np.zeros(3,'i')

  offset = min(options.z)/info['size'][2]*info['grid'][2]                                   # offset due to lower z-plane
  n = 0
  for i in xrange(Nx):
    grid[0] = round((i+0.5)*info['grid'][0]/Nx-0.5)
    for j in xrange(Ny):
      grid[1] = round((j+0.5)*info['grid'][1]/Ny-0.5)
      for k in xrange(Nz):
        grid[2] = offset + k
        grid %= info['grid']
        coordinates = (0.5+grid)*info['size']/info['grid']
        seeds[n,0:3] = coordinates/info['size']                                           # normalize coordinates to box
        seeds[n,  3] = microstructure[grid[0],grid[1],grid[2]]
#        file['croak'].write('%s\t%i\n'%(str(seeds[n,:3]),seeds[n,3]))
        if options.x: grid[0] += 1
        if options.y: grid[1] += 1
        n += 1
#      file['croak'].write('\n')
      
  newInfo['microstructures'] = len(np.unique(seeds[:,3]))

#--- report ---------------------------------------------------------------------------------------
  if (newInfo['microstructures'] != info['microstructures']):
    file['croak'].write('--> microstructures: %i\n'%newInfo['microstructures'])

#--- write header ---------------------------------------------------------------------------------
  theTable.labels_clear()
  theTable.labels_append(['x','y','z','microstructure'])
  theTable.info_clear()
  theTable.info_append(extra_header+[
    scriptID,
    "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
    "size\tx %f\ty %f\tz %f"%(info['size'][0],info['size'][1],info['size'][2],),
    "origin\tx %f\ty %f\tz %f"%(info['origin'][0],info['origin'][1],info['origin'][2],),
    "homogenization\t%i"%info['homogenization'],
    "microstructures\t%i"%(newInfo['microstructures']),
    ])

  theTable.head_write()
  theTable.output_flush()
  theTable.data = seeds
  theTable.data_writeArray('%g')
  theTable.output_flush()
  

#--- output finalization --------------------------------------------------------------------------
  if file['name'] != 'STDIN':
    theTable.close()
    os.rename(file['name']+'_tmp',os.path.splitext(file['name'])[0] + '_poked_%ix%ix%i.seeds'%(Nx,Ny,Nz))
