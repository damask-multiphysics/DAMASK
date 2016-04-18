#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
import damask
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
translate microstructure indices (shift or substitute) and/or geometry origin.

""", version=scriptID)

parser.add_option('-o', '--origin',
                  dest = 'origin',
                  type = 'float', nargs = 3, metavar = ' '.join(['float']*3),
                  help = 'offset from old to new origin of grid')
parser.add_option('-m', '--microstructure',
                  dest = 'microstructure',
                  type = 'int', metavar = 'int',
                  help = 'offset from old to new microstructure indices')
parser.add_option('-s', '--substitute',
                  dest = 'substitute',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'substitutions of microstructure indices from,to,from,to,...')

parser.set_defaults(origin = (0.0,0.0,0.0),
                    microstructure = 0,
                    substitute = [],
                   )

(options, filenames) = parser.parse_args()

sub = {}
for i in xrange(len(options.substitute)/2):                                                         # split substitution list into "from" -> "to"
  sub[int(options.substitute[i*2])] = int(options.substitute[i*2+1])

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
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

  microstructure = table.microstructure_read(info['grid'])                                          # read microstructure

# --- do work ------------------------------------------------------------------------------------

  newInfo = {
             'origin':  np.zeros(3,'d'),
             'microstructures': 0,
            }

  substituted = np.copy(microstructure)
  for k, v in sub.iteritems(): substituted[microstructure==k] = v                                   # substitute microstructure indices

  substituted += options.microstructure                                                             # shift microstructure indices

  newInfo['origin'] = info['origin'] + options.origin
  newInfo['microstructures'] = substituted.max()

# --- report ---------------------------------------------------------------------------------------

  remarks = []
  if (any(newInfo['origin']          != info['origin'])):
    remarks.append('--> origin   x y z:  %s'%(' : '.join(map(str,newInfo['origin']))))
  if (    newInfo['microstructures'] != info['microstructures']):
    remarks.append('--> microstructures: %i'%newInfo['microstructures'])
  if remarks != []: damask.util.croak(remarks)

# --- write header ---------------------------------------------------------------------------------

  table.labels_clear()
  table.info_clear()
  table.info_append([
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=info['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=info['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=newInfo['origin']),
    "homogenization\t{homog}".format(homog=info['homogenization']),
    "microstructures\t{microstructures}".format(microstructures=newInfo['microstructures']),
    extra_header
    ])
  table.head_write()

# --- write microstructure information ------------------------------------------------------------

  formatwidth = int(math.floor(math.log10(microstructure.max())+1))
  table.data = substituted.reshape((info['grid'][0],info['grid'][1]*info['grid'][2]),order='F').transpose()
  table.data_writeArray('%%%ii'%(formatwidth),delimiter = ' ')

# --- output finalization --------------------------------------------------------------------------

  table.close()                                                                                     # close ASCII table
