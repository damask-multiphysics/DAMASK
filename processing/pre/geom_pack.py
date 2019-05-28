#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [geomfile(s)]', description = """
Pack ranges to "a to b" and/or multiples to "n of x".

""", version = scriptID)

(options, filenames) = parser.parse_args()


if filenames == []: filenames = [None]

for name in filenames:
  damask.util.report(scriptName,name)

  geom = damask.Geom.from_file(StringIO(''.join(sys.stdin.read())) if name is None else name)
  microstructure = geom.get_microstructure().flatten('F')
  damask.util.croak(geom)
  geom.add_comments(scriptID + ' ' + ' '.join(sys.argv[1:]))
  
  compressType = None
  former = start = -1
  reps = 0

  for current in microstructure:
    if abs(current - former) == 1 and (start - current) == reps*(former - current):
      compressType = 'to'
      reps += 1
    elif current == former and start == former:
      compressType = 'of'
      reps += 1
    else:
      if   compressType is None:
        out = geom.get_header()
      elif compressType == '.':
        out.append('{}'.format(former))
      elif compressType == 'to':
        out.append('{} to {}'.format(start,former))
      elif compressType == 'of':
        out.append('{} of {}'.format(reps,former))

      compressType = '.'
      start = current
      reps = 1

    former = current
   
  if compressType == '.':
    out.append('{}'.format(former))
  elif compressType == 'to':
    out.append('{} to {}'.format(start,former))
  elif compressType == 'of':
    out.append('{} of {}'.format(reps,former))
 
  
  if name is None:
    sys.stdout.write('\n'.join(out)+'\n')
  else:
    with open(name,'w') as f:
      f.write('\n'.join(out)+'\n')



