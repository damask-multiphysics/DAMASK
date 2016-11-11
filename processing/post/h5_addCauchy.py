#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os
import damask
import numpy as np
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID = ' '.join([scriptName, damask.version])


def getCauchy(f, p):
    """Return Cauchy stress for given f and p"""
    # [Cauchy] = (1/det(F)) * [P].[F_transpose]
    f = f.reshape((3, 3))
    p = p.reshape((3, 3))
    return 1.0/np.linalg.det(f)*np.dot(p, f.T).reshape(9)


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
desp = "Add column(s) containing Cauchy stress based on given column(s)"
desp += "of deformation gradient and first Piola--Kirchhoff stress."
parser = OptionParser(option_class=damask.extendableOption,
                      usage='%prog options [file[s]]',
                      description=desp,
                      version=scriptID)
parser.add_option('-f', '--defgrad',
                  dest='defgrad',
                  type='string', metavar='string',
                  help='heading for deformation gradient [%default]')
parser.add_option('-p', '--stress',
                  dest='stress',
                  type='string', metavar='string',
                  help='heading for first Piola--Kirchhoff stress [%default]')

parser.set_defaults(defgrad='f',
                    stress='p')

(options, filenames) = parser.parse_args()

# ----- loop over input H5 files ----- #
for name in filenames:
    try:
        h5f = damask.H5Table(name, new_file=False)
    except:
        continue
    damask.util.report(scriptName, name)

    # ----- read in data ----- #
    f = h5f.get_data("f")
    p = h5f.get_data("p")

    # ----- calculate Cauchy stress ----- #
    cauchy = [getCauchy(f_i, p_i) for f_i, p_i in zip(f, p)]

    # ----- write to HDF5 file ----- #
    cmd_log = " ".join([scriptID, name])
    h5f.add_data('Cauchy', np.array(cauchy), cmd_log=cmd_log)
