#!/usr/bin/env python
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
compress geometry files with ranges "a to b" and/or multiples "n of x".

""", version = scriptID)

(options, filenames) = parser.parse_args()

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

# --- write header ---------------------------------------------------------------------------------

  table.labels_clear()
  table.info_clear()
  table.info_append(extra_header+[
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=info['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=info['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=info['origin']),
    "homogenization\t{homog}".format(homog=info['homogenization']),
    "microstructures\t{microstructures}".format(microstructures=info['microstructures']),
    ])
  table.head_write()
  
# --- write packed microstructure information -----------------------------------------------------

  type = ''
  former = start = -1
  reps = 0

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    items = table.data
    if len(items) > 2:
      if   items[1].lower() == 'of': items = [int(items[2])]*int(items[0])
      elif items[1].lower() == 'to': items = xrange(int(items[0]),1+int(items[2]))
      else:                          items = map(int,items)
    else:                            items = map(int,items)

    for current in items:
      if current == former+1 and start+reps == former+1:   
        type = 'to'
        reps += 1
      elif current == former and start == former:
        type = 'of'
        reps += 1
      else:
        if   type == '':
          table.data = []
        elif type == '.':
          table.data = [str(former)]
        elif type == 'to':
          table.data = ['{0} to {1}'.format(former-reps+1,former)]
        elif type == 'of':
          table.data = ['{0} of {1}'.format(reps,former)]

        outputAlive = table.data_write(delimiter = ' ')                                             # output processed line

        type = '.'
        start = current
        reps = 1

      former = current

  table.data = {
                '.' : [str(former)],
                'to': ['%i to %i'%(former-reps+1,former)],
                'of': ['%i of %i'%(reps,former)],
               }[type]

  outputAlive = table.data_write(delimiter = ' ')                                                   # output processed line

# --- output finalization --------------------------------------------------------------------------

  table.close()                                                                                     # close ASCII table
