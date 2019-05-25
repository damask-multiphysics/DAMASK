#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os
import sys
import numpy as np
from scipy import ndimage
from optparse import OptionParser
from io import StringIO
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

def taintedNeighborhood(stencil,trigger=[],size=1):

  me = stencil[stencil.shape[0]//2]
  if len(trigger) == 0:
    return np.any(stencil != me)
  if me in trigger:
    trigger = set(trigger)
    trigger.remove(me)
    trigger = list(trigger)
  return np.any(np.in1d(stencil,np.array(trigger)))

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Offset microstructure index for points which see a microstructure different from themselves
(or listed as triggers) within a given (cubic) vicinity, i.e. within the region close to a grain/phase boundary.

""", version = scriptID)

parser.add_option('-v', '--vicinity',
                  dest = 'vicinity',
                  type = 'int', metavar = 'int',
                  help = 'voxel distance checked for presence of other microstructure [%default]')
parser.add_option('-m', '--microstructureoffset',
                  dest='offset',
                  type = 'int', metavar = 'int',
                  help = 'offset (positive or negative) for tagged microstructure indices. '+
                         '"0" selects maximum microstructure index [%default]')
parser.add_option('-t', '--trigger',
                  action = 'extend', dest = 'trigger', metavar = '<int LIST>',
                  help = 'list of microstructure indices triggering a change')
parser.add_option('-n', '--nonperiodic',
                  dest = 'mode',
                  action = 'store_const', const = 'nearest',
                  help = 'assume geometry to be non-periodic')

parser.set_defaults(vicinity = 1,
                    offset   = 0,
                    trigger  = [],
                    mode     = 'wrap',
                   )

(options, filenames) = parser.parse_args()
options.trigger = np.array(options.trigger, dtype=int)


# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  if name is None:
    virt_file = StringIO(''.join(sys.stdin.read()))
    geom = damask.Geom.from_file(virt_file)
  else:
    try: # really needed? Why not simply fail if file does not exists etc.
      geom = damask.Geom.from_file(name)
    except: continue
  damask.util.report(scriptName,name)

  microstructure = geom.microstructure

  if options.offset == 0: options.offset = microstructure.max()

  microstructure = np.where(ndimage.filters.generic_filter(microstructure,
                                                           taintedNeighborhood,
                                                           size=1+2*options.vicinity,mode=options.mode,
                                                           extra_arguments=(),
                                                           extra_keywords={"trigger":options.trigger,"size":1+2*options.vicinity}),
                            microstructure + options.offset,microstructure)

  geom.microstructure = microstructure
  geom.add_comment(scriptID + ' ' + ' '.join(sys.argv[1:]))
  
  damask.util.croak('\n'.join(geom.info()))
  
  if name is None:
    sys.stdout.write(str(geom))
  else:
    geom.to_file(name)
