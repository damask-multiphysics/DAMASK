#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser
from io import StringIO

import numpy as np

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

sub = {}
for i in range(len(options.substitute)//2):                                                         # split substitution list into "from" -> "to"
  sub[int(options.substitute[i*2])] = int(options.substitute[i*2+1])


if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)
  
  if name is None:
    virt_file = StringIO(''.join(sys.stdin.read()))
    geom = damask.Geom.from_file(virt_file)
  else:
    geom = damask.Geom.from_file(name)
  damask.util.croak(geom)  
  microstructure = geom.get_microstructure()
  new            = np.copy(microstructure)
  
  for k, v in sub.items(): new[microstructure==k] = v                                               # substitute microstructure indices and shift
  
  microstructure += options.microstructure                                                          # constant shift

  for i,line in enumerate(geom.comments):
    if line.lower().strip().startswith('origin'):
      origin= np.array([float(line.split()[j]) for j in [2,4,6]])                                   # assume correct order (x,y,z)
      origin += np.array(origin)
      geom.comments[i] = 'origin x {} y {} z {}'.format(*origin)
  
  damask.util.croak(geom.update(microstructure))
  geom.add_comment(scriptID + ' ' + ' '.join(sys.argv[1:]))

  if name is None:
    sys.stdout.write(str(geom.show()))
  else:
    geom.to_file(name)
