# -*- coding: UTF-8 no BOM -*-

import os
import sys
import h5py
import numpy as np
import xml.etree.cElementTree as ET

# ---------------------------------------------------------------- #
# python 3 has no unicode object, this ensures that the code works #
# on Python 2&3                                                    #
# ---------------------------------------------------------------- #
try:
  test=isinstance('test', unicode)
except(NameError):
  unicode=str


# ------------------------------------------------------- #
# Singleton class for converting feature name to H5F path #
# ------------------------------------------------------- #
# NOTE:
#   use simple function to mimic the singleton class in
#   C++/Java
def lables_to_path(label, dsXMLPath=None):
    """ read the xml definition file and return the path."""
    if dsXMLPath is None:
        dsXMLPath = os.path.abspath(__file__).replace("h5table.py",
                                                      "DS_HDF5.xml")
    tree = ET.parse(dsXMLPath)
    dataType = tree.find('{}/type'.format(label)).text
    h5path = tree.find('{}/h5path'.format(label)).text
    return (dataType, h5path)


# ----------------------- #
# H5Table interface class #
# ----------------------- #
class H5Table(object):
    """
    DESCRIPTION
    -----------
        Interface class for manipulating data in HDF5 with DAMASK
        specialized data structure.
    PARAMETERS
    ----------
    h5f_path: str
        Absolute path the HDF5 file
    METHOD
    ------
    get_attr()
    add_attr()

    get_data()
    add_data()
    NOTE
    ----
        1. As an interface class, it uses the lazy evaluation design
        that read the data only when its absolutely necessary.
    """

    def __init__(self, h5f_path):
        """
        """
        self.h5f_path = h5f_path

    def __del__(self, feature_name=None):
        dataType, h5f_path = lables_to_path(feature_name)
        h5f = h5py.File(self.h5f_path, 'a')
        del h5f[h5f_path]

    def get_attr(self, attr_name=None):
        """
        """
        h5f = h5py.File(self.h5f_path, 'r')

    def add_attr(self, ):
        """
        """
        pass

    def get_data(self, feature_name=None):
        """ extract dataset from HDF5 table and return it in a numpy array """
        dataType, h5f_path = lables_to_path(feature_name)
        h5f = h5py.File(self.h5f_path, 'r')
        h5f_dst = h5f[h5f_path]  # get the handle for target dataset(table)
        return h5f_dst.read_direct(np.zeros(h5f_dst.shape))

    def add_data(self, feature_name, dataset=None):
        """ adding new feature into existing HDF5 file """
        dataType, h5f_path = lables_to_path(feature_name)
        if dataType is not "attr":
            h5f = h5py.File(self.h5f_path, 'a')
            h5f.create_dataset(h5f_path, data=dataset)
        else:
            raise ValueError("feature {} isn't valid".format(feature_name))
