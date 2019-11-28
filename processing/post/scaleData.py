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

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Uniformly scale column values by given factor.

""", version = scriptID)

parser.add_option('-l','--label',
                  dest = 'label',
                  action = 'extend', metavar = '<string LIST>',
                  help  ='column(s) to scale')
parser.add_option('-f','--factor',
                  dest = 'factor',
                  action = 'extend', metavar='<float LIST>',
                  help = 'factor(s) per column')

parser.set_defaults(label  = [],
                   )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

if len(options.label) != len(options.factor):
    parser.error('number of column labels and factors do not match.')

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
    for i,label in enumerate(options.label):
        table.set_array(label,table.get_array(label)*float(options.factor[i]),
                        scriptID+' '+' '.join(sys.argv[1:]))

    table.to_ASCII(sys.stdout if name is None else name)
