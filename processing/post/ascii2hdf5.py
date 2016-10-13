#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

# ------------------------------------------------------------------- #
# NOTE:                                                               #
#   1. Not all output is defined in the DS_HDF5.xml, please add new   #
#      new one to the system wide definition file                     #
#             <DAMASK_ROOT>/lib/damask/DS_HDF5.xml                    #
#      or specify your own when initializing HDF5 class               #
#   2. Somehow the point cloud structure cannot be properly handled   #
#      by Xdmf, which is a descriptive wrapper for visualizing HDF5   #
#      using Paraview. The current solution is using cell structured  #
#      HDF5 so that Xdmf can describe the data shape as a rectangular #
#      mesh rather than polyvertex.                                   #
# ------------------------------------------------------------------- #

import os
import damask
import numpy as np
from optparse import OptionParser


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# ----- helper function ----- #
def get_rectMshVectors(xyz_array, posNum):
    """take in a xyz array from rectangular mesh and figure out Vx, Vy, Vz"""
    # need some improvement, and only works for rectangular grid
    v = sorted(list(set(xyz_array[:, posNum])))
    v_interval = (v[2]+v[1])/2.0 - (v[1]+v[0])/2.0
    v_start = (v[1]+v[0])/2.0 - v_interval
    v_end = (v[-1]+v[-2])/2.0 + v_interval
    V = np.linspace(v_start, v_end, len(v)+1)
    return V


# ----- MAIN ---- #
desp_msg = "Convert DAMASK ascii table to HDF5 file"
parser = OptionParser(option_class=damask.extendableOption,
                      usage='%prog options [file[s]]',
                      description = desp_msg,
                      version = scriptID)
parser.add_option('-D', '--DefinitionFile',
                  dest = 'storage definition file',
                  type = 'string',
                  metavar = 'string',
                  help = 'definition file for H5 data storage')
parser.add_option('-p',
                  '--pos', '--position',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'label of coordinates [%default]')

parser.set_defaults(DefinitionFile='default',
                    pos='pos')

(options,filenames) = parser.parse_args()

filename = filenames[0]

if options.DefinitionFile == 'default':
    defFile = None
else:
    defFile = options.DefinitionFile

# ----- read in data using DAMASK ASCII table class ----- #
asciiTable = damask.ASCIItable(name=filename, buffered=False)
asciiTable.head_read()
asciiTable.data_readArray()
incNum = int(asciiTable.data[asciiTable.label_index('inc'), 0])
fullTable = np.copy(asciiTable.data)  # deep copy all data, just to be safe
labels = asciiTable.labels()
labels_idx = [asciiTable.label_index(label) for label in labels]
featuresDim = [labels_idx[i+1] - labels_idx[i] for i in xrange(len(labels)-1)]
featuresDim.append(fullTable.shape[1] - labels_idx[-1])

# ----- figure out size and grid ----- #
pos_idx = asciiTable.label_index('pos')
xyz_array = asciiTable.data[:, pos_idx:pos_idx+3]
Vx = get_rectMshVectors(xyz_array, 0)
Vy = get_rectMshVectors(xyz_array, 1)
Vz = get_rectMshVectors(xyz_array, 2)
# use the dimension of the rectangular grid to reshape all other data
mshGridDim = (len(Vx)-1, len(Vy)-1, len(Vz)-1)

# ----- create a new HDF5 file and save the data -----#
# Will overwrite existing HDF5 file with the same name
h5f = damask.H5Table(filename.replace(".txt", ".h5"),
                     new_file=True,
                     dsXMLFile=defFile)
# adding increment number as root level attributes
h5f.add_attr('inc', incNum)
# add the mesh grid data now
h5f.add_data("Vx", Vx)
h5f.add_data("Vy", Vy)
h5f.add_data("Vz", Vz)

# add the rest of data from table
for fi in xrange(len(labels)):
    featureName = labels[fi]
    if 'inc' in featureName: continue
    # remove trouble maker "("" and ")"
    if "(" in featureName: featureName = featureName.replace("(", "")
    if ")" in featureName: featureName = featureName.replace(")", "")
    featureIdx = labels_idx[fi]
    featureDim = featuresDim[fi]
    # grab the data hook
    dataset = fullTable[:, featureIdx:featureIdx+featureDim]
    # mapping 2D data onto a 3D rectangular mesh to get 4D data
    dataset = dataset.reshape((mshGridDim[0], mshGridDim[1], mshGridDim[2],
                               dataset.shape[1]))
    # write out data
    h5f.add_data(featureName, dataset)
