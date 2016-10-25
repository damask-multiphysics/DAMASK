# -*- coding: UTF-8 no BOM -*-

# ----------------------------------------------------------- #
# Ideally the h5py should be enough to serve as the data      #
# interface for future DAMASK, but since we are still not     #
# sure when this major shift will happen, it seems to be a    #
# good idea to provide a interface class that help user ease  #
# into using HDF5 as the new daily storage driver.            #
# ----------------------------------------------------------- #

import os
import h5py
import numpy as np
import xml.etree.cElementTree as ET

# ---------------------------------------------------------------- #
# python 3 has no unicode object, this ensures that the code works #
# on Python 2&3                                                    #
# ---------------------------------------------------------------- #
try:
    test = isinstance('test', unicode)
except(NameError):
    unicode = str


def lables_to_path(label, dsXMLPath=None):
    """Read the xml definition file and return the path."""
    if dsXMLPath is None:
        # use the default storage layout in DS_HDF5.xml
        if "h5table.pyc" in __file__:
            dsXMLPath = os.path.abspath(__file__).replace("h5table.pyc",
                                                          "DS_HDF5.xml")
        else:
            dsXMLPath = os.path.abspath(__file__).replace("h5table.py",
                                                          "DS_HDF5.xml")
    # This current implementation requires that all variables
    # stay under the root node, the nesting is defined through the
    # h5path.
    # Allow new derived data to be put under the root
    tree = ET.parse(dsXMLPath)
    try:
        dataType = tree.find('{}/type'.format(label)).text
        h5path = tree.find('{}/h5path'.format(label)).text
    except:
        dataType = "Scalar"
        h5path = "/{}".format(label)  # just put it under root
    return (dataType, h5path)


class H5Table(object):
    """
    Light weight interface class for h5py

    DESCRIPTION
    -----------
        Interface/wrapper class for manipulating data in HDF5 with DAMASK
        specialized data structure.
        -->try to maintain a minimal API design.
    PARAMETERS
    ----------
    h5f_path: str
        Absolute path the HDF5 file
    METHOD
    ------
    del_entry()  -- Force delete attributes/group/datasets (Dangerous)
    get_attr()   -- Return attributes if possible
    add_attr()   -- Add NEW attributes to dataset/group (no force overwrite)
    get_data()   -- Retrieve data in numpy.ndarray
    add_data()   -- Add dataset to H5 file
    get_cmdlog() -- Return the command used to generate the data if possible.
    NOTE
    ----
        1. As an interface class, it uses the lazy evaluation design
        that read the data only when its absolutely necessary.
        2. The command line used to generate new feature is stored with
        each dataset as dataset attribute.

    """

    def __init__(self, h5f_path, new_file=False, dsXMLFile=None):
        self.h5f_path = h5f_path
        self.dsXMLFile = dsXMLFile
        msg = 'Created by H5Talbe from DAMASK'
        mode = 'w' if new_file else 'a'
        with h5py.File(self.h5f_path, mode) as h5f:
            h5f['/'].attrs['description'] = msg

    def del_entry(self, feature_name):
        """Delete entry in HDF5 table"""
        dataType, h5f_path = lables_to_path(feature_name,
                                            dsXMLPath=self.dsXMLFile)
        with h5py.File(self.h5f_path, 'a') as h5f:
            del h5f[h5f_path]

    def get_attr(self, attr_name):
        dataType, h5f_path = lables_to_path(attr_name,
                                            dsXMLPath=self.dsXMLFile)
        with h5py.File(self.h5f_path, 'a') as h5f:
            rst_attr = h5f[h5f_path].attrs[attr_name]
        return rst_attr

    def add_attr(self, attr_name, attr_data):
        dataType, h5f_path = lables_to_path(attr_name,
                                            dsXMLPath=self.dsXMLFile)
        with h5py.File(self.h5f_path, 'a') as h5f:
            h5f[h5f_path].attrs[attr_name] = attr_data
            h5f.flush()

    def get_data(self, feature_name=None):
        """Extract dataset from HDF5 table and return it in a numpy array"""
        dataType, h5f_path = lables_to_path(feature_name,
                                            dsXMLPath=self.dsXMLFile)
        with h5py.File(self.h5f_path, 'a') as h5f:
            h5f_dst = h5f[h5f_path]  # get the handle for target dataset(table)
            rst_data = np.zeros(h5f_dst.shape)
            h5f_dst.read_direct(rst_data)
        return rst_data

    def add_data(self, feature_name, dataset, cmd_log=None):
        """Adding new feature into existing HDF5 file"""
        dataType, h5f_path = lables_to_path(feature_name,
                                            dsXMLPath=self.dsXMLFile)
        with h5py.File(self.h5f_path, 'a') as h5f:
            # NOTE:
            # --> If dataset exists, delete the old one so as to write
            #     a new one. For brand new dataset. For brand new one,
            #     record its state as fresh in the cmd log.
            try:
                del h5f[h5f_path]
                print("***deleting old {} from {}".format(feature_name,self.h5f_path))
            except:
                # if no cmd log, None will used
                cmd_log = str(cmd_log) + " [FRESH]"
            h5f.create_dataset(h5f_path, data=dataset)
            # store the cmd in log is possible
            if cmd_log is not None:
                h5f[h5f_path].attrs['log'] = str(cmd_log)
            h5f.flush()

    def get_cmdlog(self, feature_name):
        """Get cmd history used to generate the feature"""
        dataType, h5f_path = lables_to_path(feature_name,
                                            dsXMLPath=self.dsXMLFile)
        with h5py.File(self.h5f_path, 'a') as h5f:
            cmd_logs = h5f[h5f_path].attrs['log']
        return cmd_logs
