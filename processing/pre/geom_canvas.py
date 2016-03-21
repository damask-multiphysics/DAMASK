#!/usr/bin/env python
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

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Changes the (three-dimensional) canvas of a spectral geometry description.

""", version = scriptID)

parser.add_option('-g', '--grid',
                  dest = 'grid',
                  type = 'string', nargs = 3, metavar = ' '.join(['string']*3),
                  help = 'a,b,c grid of hexahedral box [auto]')
parser.add_option('-o', '--offset',
                  dest = 'offset',
                  type = 'int', nargs = 3, metavar = ' '.join(['int']*3),
                  help = 'a,b,c offset from old to new origin of grid [%default]')
parser.add_option('-f', '--fill',
                  dest = 'fill',
                  type = 'float', metavar = 'float',
                  help = '(background) canvas grain index. "0" selects maximum microstructure index + 1 [%default]')
parser.add_option('--float',
                  dest = 'real',
                  action = 'store_true',
                  help = 'input data is float [%default]')

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
  xindex = list(set(xrange(options.offset[0],options.offset[0]+newInfo['grid'][0])) & \
                                                               set(xrange(info['grid'][0])))
  yindex = list(set(xrange(options.offset[1],options.offset[1]+newInfo['grid'][1])) & \
                                                               set(xrange(info['grid'][1])))
  zindex = list(set(xrange(options.offset[2],options.offset[2]+newInfo['grid'][2])) & \
                                                               set(xrange(info['grid'][2])))
  translate_x = [i - options.offset[0] for i in xindex]
  translate_y = [i - options.offset[1] for i in yindex]
  translate_z = [i - options.offset[2] for i in zindex]
  if 0 in map(len,[xindex,yindex,zindex,translate_x,translate_y,translate_z]):
    damask.util.croak('Invaldid grid-offset comination')
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
    remarks.append('--> grid     a b c:  %s'%(' x '.join(map(str,newInfo['grid']))))
  if (any(newInfo['size']            != info['size'])):
    remarks.append('--> size     x y z:  %s'%(' x '.join(map(str,newInfo['size']))))
  if (any(newInfo['origin']          != info['origin'])):
    remarks.append('--> origin   x y z:  %s'%(' : '.join(map(str,newInfo['origin']))))
  if (    newInfo['microstructures'] != info['microstructures']):
    remarks.append('--> microstructures: %i'%newInfo['microstructures'])

  if np.any(newInfo['grid'] < 1):    errors.append('invalid new grid a b c.')
  if np.any(newInfo['size'] <= 0.0): errors.append('invalid new size x y z.')

  if remarks != []: damask.util.croak(remarks)
  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --- write header ---------------------------------------------------------------------------------

  table.info_clear()
  table.info_append([
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=newInfo['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=newInfo['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=newInfo['origin']),
    "homogenization\t{homog}".format(homog=info['homogenization']),
    "microstructures\t{microstructures}".format(microstructures=newInfo['microstructures']),
    extra_header
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
