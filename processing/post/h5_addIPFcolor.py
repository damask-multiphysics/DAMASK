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

# TODO
#   This implementation will have to iterate through the array one
#   element at a time, maybe there are some other ways to make this
#   faster.

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
desp = "Add RGB color value corresponding to TSL-OIM scheme for IPF."
parser = OptionParser(option_class=damask.extendableOption,
                      usage='%prog options [file[s]]',
                      description=desp,
                      version=scriptID)
parser.add_option('-p', '--pole',
                  dest='pole',
                  type='float', nargs=3, metavar='float float float',
                  help='lab frame direction for IPF [%default]')
msg = ', '.join(damask.Symmetry.lattices[1:])
parser.add_option('-s', '--symmetry',
                  dest='symmetry',
                  type='choice', choices=damask.Symmetry.lattices[1:],
                  metavar='string',
                  help='crystal symmetry [%default] {{{}}} '.format(msg))
parser.add_option('-e', '--eulers',
                  dest='eulers',
                  type='string', metavar='string',
                  help='Euler angles label')
parser.add_option('-d', '--degrees',
                  dest='degrees',
                  action='store_true',
                  help='Euler angles are given in degrees [%default]')
parser.add_option('-m', '--matrix',
                  dest='matrix',
                  type='string', metavar='string',
                  help='orientation matrix label')
parser.add_option('-a',
                  dest='a',
                  type='string', metavar='string',
                  help='crystal frame a vector label')
parser.add_option('-b',
                  dest='b',
                  type='string', metavar='string',
                  help='crystal frame b vector label')
parser.add_option('-c',
                  dest='c',
                  type='string', metavar='string',
                  help='crystal frame c vector label')
parser.add_option('-q', '--quaternion',
                  dest='quaternion',
                  type='string', metavar='string',
                  help='quaternion label')

parser.set_defaults(pole=(0.0, 0.0, 1.0),
                    symmetry=damask.Symmetry.lattices[-1],
                    degrees=False)

(options, filenames) = parser.parse_args()

# safe guarding to have only one orientation representation
# use dynamic typing to group a,b,c into frame
options.frame = [options.a, options.b, options.c]
input = [options.eulers is not None,
         all(options.frame),
         options.matrix is not None,
         options.quaternion is not None]

if np.sum(input) != 1:
    parser.error('needs exactly one input format.')

# select input label that was requested (active)
label_active = np.where(input)[0][0]
(label, dim, inputtype) = [(options.eulers, 3, 'eulers'),
                           (options.frame, [3, 3, 3], 'frame'),
                           (options.matrix, 9, 'matrix'),
                           (options.quaternion, 4, 'quaternion')][label_active]

# rescale degrees to radians
toRadians = math.pi/180.0 if options.degrees else 1.0

# only use normalized pole
pole = np.array(options.pole)
pole /= np.linalg.norm(pole)

# ----- Loop over input files ----- #
for name in filenames:
    try:
        h5f = damask.H5Table(name, new_file=False)
    except:
        continue
    damask.util.report(scriptName, name)

    # extract data from HDF5 file
    if inputtype == 'eulers':
        orieData = h5f.get_data(label)
    elif inputtype == 'matrix':
        orieData = h5f.get_data(label)
        orieData = orieData.reshape(orieData.shape[0], 3, 3)
    elif inputtype == 'frame':
        vctr_a = h5f.get_data(label[0])
        vctr_b = h5f.get_data(label[1])
        vctr_c = h5f.get_data(label[2])
        frame = np.column_stack((vctr_a, vctr_b, vctr_c))
        orieData = frame.reshape(frame.shape[0], 3, 3)
    elif inputtype == 'quaternion':
        orieData = h5f.get_data(label)

    # calculate the IPF color
    rgbArrays = np.zeros((orieData.shape[0], 3))
    for ci in xrange(rgbArrays.shape[0]):
        if inputtype == 'eulers':
            o = damask.Orientation(Eulers=np.array(orieData[ci, :])*toRadians,
                                   symmetry=options.symmetry).reduced()
        elif inputtype == 'matrix':
            o = damask.Orientation(matrix=orieData[ci, :, :].transpose(),
                                   symmetry=options.symmetry).reduced()
        elif inputtype == 'frame':
            o = damask.Orientation(matrix=orieData[ci, :, :],
                                   symmetry=options.symmetry).reduced()
        elif inputtype == 'quaternion':
            o = damask.Orientation(quaternion=orieData[ci, :],
                                   symmetry=options.symmetry).reduced()
        rgbArrays[ci, :] = o.IPFcolor(pole)

    # compose labels/headers for IPF color (RGB)
    labelIPF = 'IPF_{:g}{:g}{:g}_{sym}'.format(*options.pole,
                                               sym=options.symmetry.lower())

    # compose cmd history (go with dataset)
    cmd_log = scriptID + '\t' + ' '.join(sys.argv[1:])

    # write data to HDF5 file
    h5f.add_data(labelIPF, rgbArrays, cmd_log=cmd_log)
