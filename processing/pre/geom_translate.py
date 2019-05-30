#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser
from io import StringIO

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [geomfile(s)]', description = """
Translate microstructure indices (shift or substitute) and/or geometry origin.

""", version=scriptID)

parser.add_option('-o', '--origin',
                  dest = 'origin',
                  type = 'float', nargs = 3, metavar = ' '.join(['float']*3),
                  help = 'offset from old to new origin of grid')
parser.add_option('-m', '--microstructure',
                  dest = 'microstructure',
                  type = 'int', metavar = 'int',
                  help = 'offset from old to new microstructure indices (after substitution)')
parser.add_option('-s', '--substitute',
                  dest = 'substitute',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'substitutions of microstructure indices from,to,from,to,...')

parser.set_defaults(origin = (0.0,0.0,0.0),
                    microstructure = 0,
                    substitute = []
                   )

(options, filenames) = parser.parse_args()

sub = list(map(int,options.substitute))


if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)

  geom = damask.Geom.from_file(StringIO(''.join(sys.stdin.read())) if name is None else name)

  substituted = geom.get_microstructure()
  for old,new in zip(sub[0::2],sub[1::2]): substituted[substituted==old] = new                      # substitute microstructure indices
  substituted += options.microstructure                                                             # constant shift

  damask.util.croak(geom.update(substituted,origin=geom.get_origin()+options.origin))
  geom.add_comments(scriptID + ' ' + ' '.join(sys.argv[1:]))

  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
