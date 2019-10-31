#!/usr/bin/env python3

import os
from optparse import OptionParser

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Add column containing Cauchy stress based on deformation gradient and first Piola--Kirchhoff stress.

""", version = scriptID)

parser.add_option('-f','--defgrad',
                  dest = 'defgrad',
                  type = 'string', metavar = 'string',
                  help = 'heading of columns containing deformation gradient [%default]')
parser.add_option('-p','--stress',
                  dest = 'stress',
                  type = 'string', metavar = 'string',
                  help = 'heading of columns containing first Piola--Kirchhoff stress [%default]')

parser.set_defaults(defgrad = 'f',
                    stress  = 'p',
                   )

(options,filenames) = parser.parse_args()

for name in filenames:
    table = damask.Table(name)
    table.add_array('Cauchy',damask.mechanics.Cauchy(table.get_array(options.defgrad).reshape(-1,3,3),
                                                     table.get_array(options.stress).reshape(-1,3,3)).reshape(-1,9),
                    scriptID)
    table.to_ASCII()
