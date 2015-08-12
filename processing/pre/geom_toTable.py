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

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [geomfile[s]]', description = """
Produce ASCIItable of structure data from geom description

""", version = scriptID)

parser.add_option('-p','--position',
                  dest = 'position',
                  type = 'string', metavar = 'string',
                  help = 'column label for position [%default]')

parser.set_defaults(position = 'pos',
                   )

(options, filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = ['STDIN']

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              outname = os.path.splitext(name)[0]+'.txt' if name else name,
                              buffered = False, labeled = False)
  except: continue
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name else ''))

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()
  info,extra_header = table.head_getGeom()
  
  table.croak(['grid     a b c:  %s'%(' x '.join(map(str,info['grid']))),
               'size     x y z:  %s'%(' x '.join(map(str,info['size']))),
               'origin   x y z:  %s'%(' : '.join(map(str,info['origin']))),
               'homogenization:  %i'%info['homogenization'],
               'microstructures: %i'%info['microstructures'],
              ])

  errors = []
  if np.any(info['grid'] < 1):    errors.append('invalid grid a b c.')
  if np.any(info['size'] <= 0.0): errors.append('invalid size x y z.')
  if errors != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue

# --- read data ------------------------------------------------------------------------------------

  microstructure = table.microstructure_read(info['grid'])                                          # read microstructure

# ------------------------------------------ assemble header ---------------------------------------

  table.info_clear()
  table.info_append(extra_header + [scriptID + '\t' + ' '.join(sys.argv[1:])])
  table.labels_clear()
  table.labels_append(['{dim}_{label}'.format(dim = 1+i,label = options.position) for i in range(3)]+['microstructure'])
  table.head_write()
  table.output_flush()

#--- generate grid --------------------------------------------------------------------------------

  x = (0.5 + np.arange(info['grid'][0],dtype=float))/info['grid'][0]*info['size'][0]+info['origin'][0]
  y = (0.5 + np.arange(info['grid'][1],dtype=float))/info['grid'][1]*info['size'][1]+info['origin'][1]
  z = (0.5 + np.arange(info['grid'][2],dtype=float))/info['grid'][2]*info['size'][2]+info['origin'][2]

  xx = np.tile(          x,                info['grid'][1]* info['grid'][2])
  yy = np.tile(np.repeat(y,info['grid'][0]                ),info['grid'][2])
  zz =         np.repeat(z,info['grid'][0]*info['grid'][1])

  table.data = np.squeeze(np.dstack((xx,yy,zz,microstructure)))
  table.data_writeArray()

# ------------------------------------------ finalize output ---------------------------------------

  table.close()
