#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import numpy as np

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [geomfile(s)]', description = """
Offset microstructure index for points which see a microstructure different from themselves
(or listed as triggers) within a given (cubic) vicinity, i.e. within the region close to a grain/phase boundary.

""", version = scriptID)

parser.add_option('-v', '--vicinity',
                  dest = 'vicinity',
                  type = 'int', metavar = 'int',
                  help = 'voxel distance checked for presence of other microstructure [%default]')
parser.add_option('-o', '--offset',
                  dest='offset',
                  type = 'int', metavar = 'int',
                  help='offset (positive or negative) to tag microstructure indices, defaults to max microstructure index')
parser.add_option('-t', '--trigger',
                  dest = 'trigger',
                  action = 'extend', metavar = '<int LIST>',
                  help = 'list of microstructure indices triggering a change')
parser.add_option('-n', '--nonperiodic',
                  dest = 'mode',
                  action = 'store_const', const = 'nearest',
                  help = 'assume geometry to be non-periodic')

parser.set_defaults(vicinity = 1,
                    trigger  = [],
                    mode     = 'wrap',
                   )

(options, filenames) = parser.parse_args()

options.trigger = np.array(options.trigger, dtype=int)


if filenames == []: filenames = [None]

for name in filenames:
    damask.util.report(scriptName,name)

    geom = damask.Geom.from_file(StringIO(''.join(sys.stdin.read())) if name is None else name)

    damask.util.croak(geom.vicinity_offset(options.vicinity,options.offset,options.trigger,
                                           True if options.mode is 'wrap' else False))
    geom.add_comments(scriptID + ' ' + ' '.join(sys.argv[1:]))

    geom.to_file(sys.stdout if name is None else name,pack=False)
