#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os
import sys
import math
import damask
import numpy as np
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID = ' '.join([scriptName, damask.version])


# ----- Helper functions ----- #
def calcMises(what, tensor):
    """Calculate von Mises equivalent"""
    dev = tensor - np.trace(tensor)/3.0*np.eye(3)
    symdev = 0.5*(dev+dev.T)
    return math.sqrt(np.sum(symdev*symdev.T) *
                     {
                     'stress': 3.0/2.0,
                     'strain': 2.0/3.0,
                     }[what.lower()])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
desp = "Add von Mises equivalent values for symmetric part of requested"
parser = OptionParser(option_class=damask.extendableOption,
                      usage='%prog options [file[s]]',
                      description=desp,
                      version=scriptID)
parser.add_option('-e', '--strain',
                  dest='strain',
                  metavar='string',
                  help='name of dataset containing strain tensors')
parser.add_option('-s', '--stress',
                  dest='stress',
                  metavar='string',
                  help='name of dataset containing stress tensors')

parser.set_defaults(strain=None, stress=None)

(options, filenames) = parser.parse_args()

# ----- Loop over input files ----- #
for name in filenames:
    try:
        h5f = damask.H5Table(name, new_file=False)
    except:
        continue
    damask.util.report(scriptName, name)

    # TODO:
    #   Could use some refactoring here
    if options.stress is not None:
        # extract stress tensor from HDF5
        tnsr = h5f.get_data(options.stress)

        # calculate von Mises equivalent row by row
        vmStress = np.zeros(tnsr.shape[0])
        for ri in range(tnsr.shape[0]):
            stressTnsr = tnsr[ri, :].reshape(3, 3)
            vmStress[ri] = calcMises('stress', stressTnsr)

        # compose label
        label = "Mises{}".format(options.stress)

        # prepare log info
        cmd_log = scriptID + '\t' + ' '.join(sys.argv[1:])

        # write data to HDF5 file
        h5f.add_data(label, vmStress, cmd_log=cmd_log)

    if options.strain is not None:
        tnsr = h5f.get_data(options.strain)
        vmStrain = np.zeros(tnsr.shape[0])
        for ri in range(tnsr.shape[0]):
            strainTnsr = tnsr[ri, :].reshape(3, 3)
            vmStrain[ri] = calcMises('strain', strainTnsr)
        label = "Mises{}".format(options.strain)
        cmd_log = scriptID + '\t' + ' '.join(sys.argv[1:])
        h5f.add_data(label, vmStrain, cmd_log=cmd_log)
