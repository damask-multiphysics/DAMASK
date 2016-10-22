#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os
# import re
# import sys
import collections
# import math
import damask
# import numpy as np
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID = ' '.join([scriptName, damask.version])


# ----- Helper functions ----- #
def listify(x):
    return x if isinstance(x, collections.Iterable) else [x]


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
usageEx = """
usage_in_details:
    Column labels are tagged by '#label#' in formulas.
    Use ';' for ',' in functions. Numpy is available as 'np'.
    Special variables: #_row_# -- row index

    Examples:
    (1) magnitude of vector -- "np.linalg.norm(#vec#)"
    (2) rounded root of row number -- "round(math.sqrt(#_row_#);3)"
"""
desp = "Add or alter column(s) with derived values according to "
desp += "user-defined arithmetic operation between column(s)."

parser = OptionParser(option_class=damask.extendableOption,
                      usage='%prog options [file[s]]' + usageEx,
                      description=desp,
                      version=scriptID)
parser.add_option('-l', '--label',
                  dest='labels',
                  action='extend', metavar='<string LIST>',
                  help='(list of) new column labels')
parser.add_option('-f', '--formula',
                  dest='formulas',
                  action='extend', metavar='<string LIST>',
                  help='(list of) formulas corresponding to labels')
parser.add_option('-c', '--condition',
                  dest='condition', metavar='string',
                  help='condition to filter rows')

parser.set_defaults(condition=None)

(options, filenames) = parser.parse_args()

# ----- parse formulas ----- #
for i in xrange(len(options.formulas)):
    options.formulas[i] = options.formulas[i].replace(';', ',')

# ----- loop over input files ----- #
for name in filenames:
    try:
        h5f = damask.H5Table(name, new_file=False)
    except:
        print "!!!Cannot process {}".format(name)
        continue
    damask.util.report(scriptName, name)

# Note:
# --> not immediately needed, come back later
