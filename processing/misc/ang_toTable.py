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

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [angfile[s]]', description = """
Convert TSL/EDAX *.ang file to ASCIItable

""", version = scriptID)

(options, filenames) = parser.parse_args()
if filenames == []: filenames = [None]

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.from_ang(StringIO(''.join(sys.stdin.read())) if name is None else name)
    table.to_ASCII(sys.stdout if name is None else os.path.splitext(name)[0]+'.txt')
