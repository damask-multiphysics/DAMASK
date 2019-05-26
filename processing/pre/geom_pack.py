#!/usr/bin/env python3

import os
import sys
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

  if name is None:
    virt_file = StringIO(''.join(sys.stdin.read()))
    geom = damask.Geom.from_file(virt_file)
  else:
    geom = damask.Geom.from_file(name)
  damask.util.croak(geom)
  microstructure = geom.get_microstructure().flatten('F')

  
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
      if   compressType == None:
        out = []
      elif compressType == '.':
        out.append('{}\n'.format(former))
      elif compressType == 'to':
        out.append('{} to {}\n'.format(start,former))
      elif compressType == 'of':
        out.append('{} of {}\n'.format(reps,former))

      compressType = '.'
      start = current
      reps = 1

    former = current
   
  if compressType == '.':
    out.append('{}\n'.format(former))
  elif compressType == 'to':
    out.append('{} to {}\n'.format(start,former))
  elif compressType == 'of':
    out.append('{} of {}\n'.format(reps,former))
 
  geom.add_comment(scriptID + ' ' + ' '.join(sys.argv[1:]))
  comments = geom.get_comments()
  with open(name,'w') as f:
    f.write('{} header\n'.format(3+len(comments)))
    f.writelines(["{}\n".format(comment) for comment in comments])
    f.write('grid a {} b {} c {}\n'.format(*geom.get_grid()))
    f.write('size x {} y {} z {}\n'.format(*geom.get_size()))
    f.write('homogenization {}\n'.format(geom.get_homogenization()))
    f.writelines(out)
