#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog option(s) [geomfile(s)]', description = """
Changes the (three-dimensional) canvas of a spectral geometry description.
Grid can be given as absolute or relative values, e.g. 16 16 16 or 2x 0.5x 32.

""", version = scriptID)

parser.add_option('-g',
                  '--grid',
                  dest = 'grid',
                  type = 'string', nargs = 3, metavar = ' '.join(['string']*3),
                  help = 'a,b,c grid of hexahedral box. [auto]')
parser.add_option('-o',
                  '--offset',
                  dest = 'offset',
                  type = 'int', nargs = 3, metavar = ' '.join(['int']*3),
                  help = 'a,b,c offset from old to new origin of grid [%default]')
parser.add_option('-f',
                  '--fill',
                  dest = 'fill',
                  type = 'float', metavar = 'float',
                  help = '(background) canvas grain index. "0" selects maximum microstructure index + 1 [%default]')
parser.add_option('--float',
                  dest = 'real',
                  action = 'store_true',
                  help = 'use float input')

parser.set_defaults(grid = ['0','0','0'],
                    offset = (0,0,0),
                    fill = 0,
                    real = False,
                   )

(options, filenames) = parser.parse_args()

datatype = 'f' if options.real else 'i'


# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try: table = damask.ASCIItable(name = name,
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

  microstructure = table.microstructure_read(info['grid'],datatype).reshape(info['grid'],order='F')          # read microstructure

# --- do work ------------------------------------------------------------------------------------
 
  newInfo = {
             'grid':    np.zeros(3,'i'),
             'origin':  np.zeros(3,'d'),
             'microstructures': 0,
            }

  newInfo['grid'] = np.array([int(o*float(n.translate(None,'xX'))) if n[-1].lower() == 'x'\
                                                                   else int(n) for o,n in zip(info['grid'],options.grid)],'i')
  newInfo['grid'] = np.where(newInfo['grid'] > 0, newInfo['grid'],info['grid'])

  microstructure_cropped = np.zeros(newInfo['grid'],datatype)
  microstructure_cropped.fill(options.fill if options.real or options.fill > 0 else microstructure.max()+1)
  xindex = list(set(range(options.offset[0],options.offset[0]+newInfo['grid'][0])) & \
                                                               set(range(info['grid'][0])))
  yindex = list(set(range(options.offset[1],options.offset[1]+newInfo['grid'][1])) & \
                                                               set(range(info['grid'][1])))
  zindex = list(set(range(options.offset[2],options.offset[2]+newInfo['grid'][2])) & \
                                                               set(range(info['grid'][2])))
  translate_x = [i - options.offset[0] for i in xindex]
  translate_y = [i - options.offset[1] for i in yindex]
  translate_z = [i - options.offset[2] for i in zindex]
  if 0 in map(len,[xindex,yindex,zindex,translate_x,translate_y,translate_z]):
    damask.util.croak('invaldid grid-offset combination.')
    table.close(dismiss = True)
    continue
  microstructure_cropped[min(translate_x):(max(translate_x)+1),\
                         min(translate_y):(max(translate_y)+1),\
                         min(translate_z):(max(translate_z)+1)] \
        = microstructure[min(xindex):(max(xindex)+1),\
                         min(yindex):(max(yindex)+1),\
                         min(zindex):(max(zindex)+1)]

  newInfo['size']   = info['size']/info['grid']*newInfo['grid']
  newInfo['origin'] = info['origin']+info['size']/info['grid']*options.offset
  newInfo['microstructures'] = microstructure_cropped.max()

# --- report ---------------------------------------------------------------------------------------

  remarks = []
  errors = []

  if (any(newInfo['grid']            != info['grid'])):
    remarks.append('--> grid     a b c:  {}'.format(' x '.join(map(str,newInfo['grid']))))
  if (any(newInfo['size']            != info['size'])):
    remarks.append('--> size     x y z:  {}'.format(' x '.join(map(str,newInfo['size']))))
  if (any(newInfo['origin']          != info['origin'])):
    remarks.append('--> origin   x y z:  {}'.format(' : '.join(map(str,newInfo['origin']))))
  if (    newInfo['microstructures'] != info['microstructures']):
    remarks.append('--> microstructures: {}'.format(newInfo['microstructures']))

  if np.any(newInfo['grid'] < 1):    errors.append('invalid new grid a b c.')
  if np.any(newInfo['size'] <= 0.0): errors.append('invalid new size x y z.')

  if remarks != []: damask.util.croak(remarks)
  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --- write header ---------------------------------------------------------------------------------

  table.info_clear()
  table.info_append(extra_header+[
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {}\tb {}\tc {}".format(*newInfo['grid']),
    "size\tx {}\ty {}\tz {}".format(*newInfo['size']),
    "origin\tx {}\ty {}\tz {}".format(*newInfo['origin']),
    "homogenization\t{}".format(info['homogenization']),
    "microstructures\t{}".format(newInfo['microstructures']),
    ])
  table.labels_clear()
  table.head_write()
  table.output_flush()

# --- write microstructure information ------------------------------------------------------------

  format = '%g' if options.real else '%{}i'.format(int(math.floor(math.log10(microstructure_cropped.max())+1)))
  table.data = microstructure_cropped.reshape((newInfo['grid'][0],newInfo['grid'][1]*newInfo['grid'][2]),order='F').transpose()
  table.data_writeArray(format,delimiter=' ')
    
# --- output finalization --------------------------------------------------------------------------

  table.close()                                                                                     # close ASCII table
