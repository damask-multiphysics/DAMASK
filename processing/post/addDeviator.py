#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(2)]', description = """
Add column(s) containing deviator of requested tensor column(s).

""", version = scriptID)

parser.add_option('-t','--tensor',
                  dest = 'tensor',
                  action = 'extend', metavar='<string LIST>',
                  help = 'heading of columns containing tensor field values')
parser.add_option('-s','--spherical',
                  dest = 'spherical',
                  action = 'store_true',
                  help = 'report spherical part of tensor (hydrostatic component, pressure)')

(options,filenames) = parser.parse_args()

if options.tensor is None:
  parser.error('no data column specified...')

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
    damask.util.report(scriptName,name)
    
    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
    for tensor in options.tensor:
         table.add_array('dev({})'.format(tensor),
                         damask.mechanics.deviatoric_part(table.get_array(tensor).reshape(-1,3,3)).reshape((-1,9)),
                         scriptID)
         if options.spherical:
             table.add_array('sph({})'.format(tensor),
                             damask.mechanics.spherical_part(table.get_array(tensor).reshape(-1,3,3)),
                             scriptID)
    table.to_ASCII(sys.stdout if name is None else name)
