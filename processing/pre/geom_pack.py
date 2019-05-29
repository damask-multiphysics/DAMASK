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
  damask.util.croak(geom)
  geom.add_comments(scriptID + ' ' + ' '.join(sys.argv[1:]))
  
  compressType = None
  former = start = -1
  reps = 0

  if name is None:
    f = sys.stdout
  else:
    f= open(name,'w')

  for current in geom.microstructure.flatten('F'):
    if abs(current - former) == 1 and (start - current) == reps*(former - current):
      compressType = 'to'
      reps += 1
    elif current == former and start == former:
      compressType = 'of'
      reps += 1
    else:
      if   compressType is None:
        f.write('\n'.join(geom.get_header())+'\n')
      elif compressType == '.':
        f.write('{}\n'.format(former))
      elif compressType == 'to':
        f.write('{} to {}\n'.format(start,former))
      elif compressType == 'of':
        f.write('{} of {}\n'.format(reps,former))

      compressType = '.'
      start = current
      reps = 1

    former = current
   
  if compressType == '.':
    f.write('{}\n'.format(former))
  elif compressType == 'to':
    f.write('{} to {}\n'.format(start,former))
  elif compressType == 'of':
    f.write('{} of {}\n'.format(reps,former))
