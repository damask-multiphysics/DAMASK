#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os
import sys
import damask
import numpy as np
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID = ' '.join([scriptName, damask.version])


# ----- Helper functions ----- #
def operator(stretch, strain, eigenvalues):
    # Albrecht Bertram: Elasticity and Plasticity of Large Deformations
    # An Introduction (3rd Edition, 2012), p. 102
    return {'V#ln': np.log(eigenvalues),
            'U#ln': np.log(eigenvalues),
            'V#Biot': (np.ones(3, 'd') - 1.0/eigenvalues),
            'U#Biot': (eigenvalues - np.ones(3, 'd')),
            'V#Green': (np.ones(3, 'd') - 1.0/eigenvalues/eigenvalues)*0.5,
            'U#Green': (eigenvalues*eigenvalues - np.ones(3, 'd'))*0.5,
            }[stretch+'#'+strain]


def calcEPS(defgrads, stretchType, strainType):
    """calculate specific type of strain tensor"""
    eps = np.zeros(defgrads.shape)  # initialize container

    # TODO:
    # this loop can use some performance boost
    # (multi-threading?)
    for ri in xrange(defgrads.shape[0]):
        f = defgrads[ri, :, :].reshape(3, 3)
        U, S, Vh = np.lingalg.svd(f)
        R = np.dot(U, Vh)  # rotation of polar decomposition
        if stretchType == 'U':
            stretch = np.dot(np.linalg.inv(R), f)  # F = RU
        elif stretchType == 'V':
            stretch = np.dot(f, np.linalg.inv(R))  # F = VR

        # kill nasty noisy data
        stretch = np.where(abs(stretch) < 1e-12, 0, stretch)

        (D, V) = np.linalg.eig(stretch)
        # flip principal component with negative Eigen values
        neg = np.where(D < 0.0)
        D[neg] *= -1.
        V[:, neg] *= -1.

        # check each vector for orthogonality
        # --> brutal force enforcing orthogonal base
        #     and re-normalize
        for i, eigval in enumerate(D):
            if np.dot(V[:, i], V[:, (i+1) % 3]) != 0.0:
                V[:, (i+1) % 3] = np.cross(V[:, (i+2) % 3], V[:, i])
                V[:, (i+1) % 3] /= np.sqrt(np.dot(V[:, (i+1) % 3],
                                                  V[:, (i+1) % 3].conj()))

        # calculate requested version of strain tensor
        d = operator(stretchType, strainType, D)
        eps[ri] = (np.dot(V, np.dot(np.diag(d), V.T)).real).reshape(9)

    return eps

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
desp = "Add column(s) containing given strains based on given stretches"
parser = OptionParser(option_class=damask.extendableOption,
                      usage='%prog options [file[s]]',
                      description=desp,
                      version=scriptID)
msg = 'material strains based on right Cauchy-Green deformation, i.e., C and U'
parser.add_option('-u', '--right',
                  dest='right',
                  action='store_true',
                  help=msg)
msg = 'spatial strains based on left Cauchy--Green deformation, i.e., B and V'
parser.add_option('-v', '--left',
                  dest='left',
                  action='store_true',
                  help=msg)
parser.add_option('-0', '--logarithmic',
                  dest='logarithmic',
                  action='store_true',
                  help='calculate logarithmic strain tensor')
parser.add_option('-1', '--biot',
                  dest='biot',
                  action='store_true',
                  help='calculate biot strain tensor')
parser.add_option('-2', '--green',
                  dest='green',
                  action='store_true',
                  help='calculate green strain tensor')
msg = 'heading(s) of columns containing deformation tensor values'
parser.add_option('-f', '--defgrad',
                  dest='defgrad',
                  action='extend',
                  metavar='<string LIST>',
                  help=msg)

parser.set_defaults(right=False, left=False,
                    logarithmic=False, biot=False, green=False,
                    defgrad=['f'])

(options, filenames) = parser.parse_args()

stretches = []
strains = []

if options.right:
    stretches.append('U')
if options.left:
    stretches.append('V')

if options.logarithmic:
    strains.append('ln')
if options.biot:
    strains.append('Biot')
if options.green:
    strains.append('Green')

if options.defgrad is None:
    parser.error('no data column specified.')

# ----- Loop over input files ----- #
for name in filenames:
    try:
        h5f = damask.H5Table(name, new_file=False)
    except:
        continue
    damask.util.report(scriptName, name)

    # extract defgrads from HDF5 storage
    F = h5f.get_data(options.defgrads)

    # allow calculate multiple types of strain within the
    # same cmd call
    for stretchType in stretches:
        for strainType in strains:
            # calculate strain tensor for this type
            eps = calcEPS(F, stretchType, strainType)

            # compose labels/headers for this strain tensor
            labelsStrain = strainType + stretchType

            # prepare log info
            cmd_log = scriptID + '\t' + ' '.join(sys.argv[1:])

            # write data to HDF5 file
            h5f.add_data(labelsStrain, eps, cmd_log=cmd_log)
