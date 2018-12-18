#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os,sys
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Create seed file taking microstructure indices from given geom file.
Indices can be black-listed or white-listed.

""", version = scriptID)

parser.add_option('-w',
                  '--white',
                  action = 'extend', metavar = '<int LIST>',
                  dest   = 'whitelist',
                  help   = 'whitelist of grain IDs')
parser.add_option('-b',
                  '--black',
                  action = 'extend', metavar = '<int LIST>',
                  dest   = 'blacklist',
                  help   = 'blacklist of grain IDs')
parser.add_option('-p',
                  '--pos', '--seedposition',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'label of coordinates [%default]')

parser.set_defaults(whitelist = [],
                    blacklist = [],
                    pos       = 'pos',
                   )

(options,filenames) = parser.parse_args()

options.whitelist = list(map(int,options.whitelist))
options.blacklist = list(map(int,options.blacklist))

# --- loop over output files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    outname = os.path.splitext(name)[0]+'.seeds' if name else name,
                                    buffered = False,
                                    labeled = False)
  except: continue
  damask.util.report(scriptName,name)

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()
  info,extra_header = table.head_getGeom()
  damask.util.report_geom(info)

  errors = []
  if np.any(info['grid'] < 1):    errors.append('invalid grid a b c.')
  if np.any(info['size'] <= 0.0): errors.append('invalid size x y z.')
  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --- read data ------------------------------------------------------------------------------------

  microstructure = table.microstructure_read(info['grid'])                                          # read (linear) microstructure

# --- generate grid --------------------------------------------------------------------------------

  x = (0.5 + np.arange(info['grid'][0],dtype=float))/info['grid'][0]*info['size'][0]+info['origin'][0]
  y = (0.5 + np.arange(info['grid'][1],dtype=float))/info['grid'][1]*info['size'][1]+info['origin'][1]
  z = (0.5 + np.arange(info['grid'][2],dtype=float))/info['grid'][2]*info['size'][2]+info['origin'][2]

  xx = np.tile(          x,                info['grid'][1]* info['grid'][2])
  yy = np.tile(np.repeat(y,info['grid'][0]                ),info['grid'][2])
  zz =         np.repeat(z,info['grid'][0]*info['grid'][1])

  mask = np.logical_and(np.in1d(microstructure,options.whitelist,invert=False) if options.whitelist != []
                                                                               else np.full_like(microstructure,True,dtype=bool), 
                        np.in1d(microstructure,options.blacklist,invert=True ) if options.blacklist != []
                                                                               else np.full_like(microstructure,True,dtype=bool))

# ------------------------------------------ assemble header ---------------------------------------

  table.info_clear()
  table.info_append(extra_header+[
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {}\tb {}\tc {}".format(*info['grid']),
    "size\tx {}\ty {}\tz {}".format(*info['size']),
    "origin\tx {}\ty {}\tz {}".format(*info['origin']),
    "homogenization\t{}".format(info['homogenization']),
    "microstructures\t{}".format(info['microstructures']),
    ])
  table.labels_clear()
  table.labels_append(['{dim}_{label}'.format(dim = 1+i,label = options.pos) for i in range(3)]+['microstructure'])
  table.head_write()
  table.output_flush()
  
# --- write seeds information ------------------------------------------------------------

  table.data = np.squeeze(np.dstack((xx,yy,zz,microstructure)))[mask]
  table.data_writeArray()

# ------------------------------------------ finalize output ---------------------------------------

  table.close()
