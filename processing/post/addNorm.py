#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser

import numpy as np

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# definition of element-wise p-norms for matrices
# ToDo: better use numpy.linalg.norm

def norm(which,object):

  if which == 'Abs':                                                                                # p = 1
    return sum(map(abs, object))
  elif which == 'Frobenius':                                                                        # p = 2
    return np.sqrt(sum([x*x for x in object]))
  elif which == 'Max':                                                                              # p = inf
    return max(map(abs, object))
  else:
    return -1


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Add column(s) containing norm of requested column(s) being either vectors or tensors.

""", version = scriptID)

normChoices = ['abs','frobenius','max']
parser.add_option('-n','--norm',
                  dest = 'norm',
                  type = 'choice', choices = normChoices, metavar='string',
                  help = 'type of element-wise p-norm [frobenius] {%s}'%(','.join(map(str,normChoices))))
parser.add_option('-l','--label',
                  dest = 'labels',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'heading of column(s) to calculate norm of')

parser.set_defaults(norm = 'frobenius',
                   )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

if options.norm.lower() not in normChoices:
    parser.error('invalid norm ({}) specified.'.format(options.norm))
if options.labels is None:
    parser.error('no data column specified.')

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
    for label in options.labels:
        data = table.get(label)
        data_norm = np.empty((data.shape[0],1))
        for i,d in enumerate(data):
            data_norm[i] = norm(options.norm.capitalize(),d)

        table.add('norm{}({})'.format(options.norm.capitalize(),label),
                  data_norm,
                  scriptID+' '+' '.join(sys.argv[1:]))

    table.to_ASCII(sys.stdout if name is None else name)
