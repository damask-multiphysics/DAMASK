#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,math,sys
import numpy as np
import damask
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

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
parser.add_option('-b', '--box',
                  dest = 'box',
                  type = 'float', nargs = 6, metavar = ' '.join(['float']*6),
                  help = 'bounding box as fraction in x, y, and z directions')
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
                    box = [0.0,1.0,0.0,1.0,0.0,1.0],
                    N = 16,
                    position = 'pos',
                   )

(options,filenames) = parser.parse_args()

options.box = np.array(options.box).reshape(3,2)

# --- loop over output files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              outname = os.path.splitext(name)[-2]+'_poked_{}.seeds'.format(options.N) if name else name,
                              buffered = False, labeled = False)
  except: continue
  damask.util.report(scriptName,name)

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
  offset = (np.amin(options.box, axis=1)*info['grid']/info['size']).astype(int)
  box = np.amax(options.box, axis=1) - np.amin(options.box, axis=1)

  Nx = int(options.N/math.sqrt(options.N*info['size'][1]*box[1]/info['size'][0]/box[0]))
  Ny = int(options.N/math.sqrt(options.N*info['size'][0]*box[0]/info['size'][1]/box[1]))
  Nz = int(box[2]*info['grid'][2])

  damask.util.croak('poking {} x {} x {} in box {} {} {}...'.format(Nx,Ny,Nz,*box))

  seeds = np.zeros((Nx*Ny*Nz,4),'d')
  grid = np.zeros(3,'i')

  n = 0
  for i in xrange(Nx):
    for j in xrange(Ny):
      grid[0] = round((i+0.5)*box[0]*info['grid'][0]/Nx-0.5)+offset[0]
      grid[1] = round((j+0.5)*box[1]*info['grid'][1]/Ny-0.5)+offset[1]
      damask.util.croak('x,y coord on surface: {},{}...'.format(*grid[:2]))
      for k in xrange(Nz):
        grid[2] = k + offset[2]
        grid %= info['grid']
        seeds[n,0:3] = (0.5+grid)/info['grid']                                                     # normalize coordinates to box
        seeds[n,  3] = microstructure[grid[0],grid[1],grid[2]]
        if options.x: grid[0] += 1
        if options.y: grid[1] += 1
        n += 1
      
  newInfo['microstructures'] = len(np.unique(seeds[:,3]))

# --- report ---------------------------------------------------------------------------------------
  if (newInfo['microstructures'] != info['microstructures']): 
    damask.util.croak('--> microstructures: %i'%newInfo['microstructures'])


# ------------------------------------------ assemble header ---------------------------------------
  table.info_clear()
  table.info_append(extra_header+[
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "poking\ta {}\tb {}\tc {}".format(Nx,Ny,Nz),
    "grid\ta {}\tb {}\tc {}".format(*info['grid']),
    "size\tx {}\ty {}\tz {}".format(*info['size']),
    "origin\tx {}\ty {}\tz {}".format(*info['origin']),
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
