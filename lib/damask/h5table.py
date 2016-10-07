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
#   1. use simple diction to mimic the singleton class in
#   C++/Java
#   2. data structure
#   {<FEATURE_NAME>: (<TYPE>, <PATH>), ...}
#       <FEATURE_NAME> ==> the name of the feature, e.g. f, ...
#       <TYPE> ==> attribute, dataset
#       <PATH> ==> path_in_HDF5_file, e.g. '/f'
def lables_to_path(label, ns="https://damask.mpie.de"):
    """ read the xml definition file and return the path."""
    dsXMLPath = os.path.abspath(__file__).replace("h5table.py",
                                                  "DS_HDF5.xml")
    tree = ET.parse(dsXMLPath)
    root = tree.getroot()
    path = 'test'
    return path


# ----------------------- #
# H5Table interface class #
# ----------------------- #
class H5Table(object):
    """
    DESCRIPTION
    -----------
        Container class for interfacing with HDF5 file format. This
        is a bare-bone interface class that provide simplified getter
        and setter for retrieving and appending data and attributes.
    PARAMETERS
    ----------
    h5f_path
    METHOD
    ------
    get_attr()
    add_attr()
    del_attr()

    get_data()
    add_data()
    del_data()
    NOTE
    ----
        1. As an interface class, it uses the lazy evaluation design
        that read the data only when its absolutely necessary.
        2.
    """

    def __init__(self, h5f_path):
        """
        """
        self.h5f_path = h5f_path

    def get_attr(self, feature_name=None):
        """
        """
        h5f = h5py.File(self.h5f_path, 'r')
        pass

    def add_attr(self, ):
        """
        """
        pass

    def del_attr(self, ):
        """
        """
        pass

    def get_data(self, ):
        """
        """
        pass

    def add_data(self, ):
        """
        """
        pass

    def del_data(self, ):
        """
        """
        pass