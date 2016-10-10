#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

# ------------------------------------------------------------------- #
# NOTE:                                                               #
#   1. Not all output is defined in the DS_HDF5.xml, please add new   #
#      new one to the system wide definition file                     #
#             <DAMASK_ROOT>/lib/damask/DS_HDF5.xml                    #
#      or specify your own when initializing HDF5 class               #
# ------------------------------------------------------------------- #

import os
import damask
import numpy as np
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

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

parser.set_defaults(DefinitionFile='default')

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

# ----- create a new HDF5 file and save the data -----#
# Will overwrite existing HDF5 file with the same name

h5f = damask.H5Table(filename.replace(".txt", ".h5"),
                     new_file=True,
                     dsXMLFile=defFile)
# adding increment number as root level attributes
h5f.add_attr('inc', incNum)
for fi in xrange(len(labels)):
    featureName = labels[fi]
    if featureName == 'inc': continue
    # remove trouble maker "("" and ")"
    if "(" in featureName: featureName = featureName.replace("(", "")
    if ")" in featureName: featureName = featureName.replace(")", "")
    featureIdx = labels_idx[fi]
    featureDim = featuresDim[fi]
    # grab the data hook
    dataset = fullTable[:, featureIdx:featureIdx+featureDim]
    # reshape tensor to 3x3 --> just make life slightly easier
    if dataset.shape[1] == 9:
        dataset = dataset.reshape((dataset.shape[0], 3, 3))
    # write out data
    h5f.add_data(featureName, dataset)
