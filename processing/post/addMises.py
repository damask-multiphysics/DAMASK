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
Add vonMises equivalent values for symmetric part of requested strains and/or stresses.

""", version = scriptID)

parser.add_option('-e','--strain',
                  dest = 'strain',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'heading(s) of columns containing strain tensors')
parser.add_option('-s','--stress',
                  dest = 'stress',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'heading(s) of columns containing stress tensors')

parser.set_defaults(strain = [],
                    stress = [],
                   )
(options,filenames) = parser.parse_args()

if options.stress is [] and options.strain is []:
    parser.error('no data column specified...')

if filenames == []: filenames = [None]

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)
    for strain in options.strain:
        table.add_array('Mises({})'.format(strain),
                        damask.mechanics.Mises_strain(damask.mechanics.symmetric(table.get_array(strain).reshape(-1,3,3))),
                        scriptID+' '+' '.join(sys.argv[1:]))
    for stress in options.stress:
        table.add_array('Mises({})'.format(stress),
                        damask.mechanics.Mises_stress(damask.mechanics.symmetric(table.get_array(stress).reshape(-1,3,3))),
                        scriptID+' '+' '.join(sys.argv[1:]))
    
    table.to_ASCII(sys.stdout if name is None else name)
