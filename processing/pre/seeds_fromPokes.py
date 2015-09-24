#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,math,sys
import numpy as np
import damask
from optparse import OptionParser

scriptID = '$Id$'
scriptName = scriptID.split()[1]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Create seeds file by poking at 45 degree through given geom file.
Mimics APS Beamline 34-ID-E DAXM poking.

""", version = scriptID)

parser.add_option('-N', '--points',
                  dest = 'N',
                  type = 'int', metavar = 'int',
                  help = 'number of poking locations [%default]')
parser.add_option('-z', '--planes',
                  dest = 'z',
                  type = 'float', nargs = 2, metavar='float float',
                  help = 'top and bottom z plane')
parser.add_option('-x',
                  action = 'store_true',
                  dest   = 'x', 
                  help   = 'poke 45 deg along x')
parser.add_option('-y',
                  action = 'store_true',
                  dest   = 'y',
                  help   = 'poke 45 deg along y')
parser.add_option('-p','--position',
                  dest = 'position',
                  type = 'string', metavar = 'string',
                  help = 'column label for coordinates [%default]')

parser.set_defaults(x = False,
                    y = False,
                    N = 16,
                    position = 'pos',
                   )

(options,filenames) = parser.parse_args()

if options.z == None: parser.error('please specify top and bottom z plane...')
# --- loop over output files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              outname = os.path.splitext(name)[-2]+'_poked_{}.seeds'.format(options.N) if name else name,
                              buffered = False, labeled = False)
  except: continue
  damask.util.report(scriptName,name)

  print os.path.splitext(name)[-1]+'_poked_{}.seeds'.format(options.N)
# --- interpret header ----------------------------------------------------------------------------

  table.head_read()
  info,extra_header = table.head_getGeom()
  
  damask.util.croak(['grid     a b c:  %s'%(' x '.join(map(str,info['grid']))),
               'size     x y z:  %s'%(' x '.join(map(str,info['size']))),
               'origin   x y z:  %s'%(' : '.join(map(str,info['origin']))),
               'homogenization:  %i'%info['homogenization'],
               'microstructures: %i'%info['microstructures'],
              ])

  errors = []
  if np.any(info['grid'] < 1):    errors.append('invalid grid a b c.')
  if np.any(info['size'] <= 0.0): errors.append('invalid size x y z.')
  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --- read data ------------------------------------------------------------------------------------

  microstructure = table.microstructure_read(info['grid']).reshape(info['grid'],order='F')          # read microstructure

# --- do work ------------------------------------------------------------------------------------
 
  newInfo = {
             'microstructures': 0,
            }

  Nx = int(options.N/math.sqrt(options.N*info['size'][1]/info['size'][0]))
  Ny = int(options.N/math.sqrt(options.N*info['size'][0]/info['size'][1]))
  Nz = int((max(options.z)-min(options.z))/info['size'][2]*info['grid'][2])

  damask.util.croak('poking {} x {} x {}...'.format(Nx,Ny,Nz))

  seeds = np.zeros((Nx*Ny*Nz,4),'d')
  grid = np.zeros(3,'i')

  offset = min(options.z)/info['size'][2]*info['grid'][2]                                           # offset due to lower z-plane
  n = 0
  for i in xrange(Nx):
    grid[0] = round((i+0.5)*info['grid'][0]/Nx-0.5)
    for j in xrange(Ny):
      grid[1] = round((j+0.5)*info['grid'][1]/Ny-0.5)
      for k in xrange(Nz):
        grid[2] = offset + k
        grid %= info['grid']
        coordinates = (0.5+grid)*info['size']/info['grid']
        seeds[n,0:3] = coordinates/info['size']                                                     # normalize coordinates to box
        seeds[n,  3] = microstructure[grid[0],grid[1],grid[2]]
        if options.x: grid[0] += 1
        if options.y: grid[1] += 1
        n += 1
      
  newInfo['microstructures'] = len(np.unique(seeds[:,3]))

# --- report ---------------------------------------------------------------------------------------

  remarks = []
  if (    newInfo['microstructures'] != info['microstructures']): remarks.append('--> microstructures: %i'%newInfo['microstructures'])
  if remarks != []: damask.util.croak(remarks)

# ------------------------------------------ assemble header ---------------------------------------
  table.info_clear()
  table.info_append(extra_header+[
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "poking\ta {}\tb {}\tc {}".format(Nx,Ny,Nz),
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=info['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=info['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=info['origin']),
    "homogenization\t{}".format(info['homogenization']),
    "microstructures\t{}".format(newInfo['microstructures']),
    ])
  table.labels_clear()
  table.labels_append(['{dim}_{label}'.format(dim = 1+i,label = options.position) for i in range(3)]+['microstructure'])
  table.head_write()
  table.output_flush()
  
# --- write seeds information ------------------------------------------------------------

  table.data = seeds
  table.data_writeArray()
    
# --- output finalization --------------------------------------------------------------------------

  table.close()                                                                                     # close ASCII table
