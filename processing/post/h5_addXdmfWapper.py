#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

# ------------------------------------------------------------------- #
# NOTE:                                                               #
#   1. Current Xdmf rendering in Paraview has some memory issue where #
#      large number of polyvertices will cause segmentation fault. By #
#      default, paraview output a cell based xdmf description, which  #
#      is working for small and medium mesh (<10,000) points. Hence a #
#      rectangular mesh is used as the de facto Geometry description  #
#      here.                                                          #
#   2. Due to the unstable state Xdmf, it is safer to use port data   #
#      to VTR rather than using xdmf as interpretive layer for data   #
#      visualization.                                                 #
# ------------------------------------------------------------------- #


import damask
import h5py
import xml.etree.cElementTree as ET
from xml.dom import minidom
from damask.h5table import lables_to_path

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# ----- HELPER FUNCTIONS -----#
def addTopLvlCmt(xmlstr, topLevelCmt):
    """add top level comment to string from ET"""
    # a quick hack to add the top level comment to XML file
    # --> somehow Elementtree does not provide this functionality
    # --> by default
    strList = xmlstr.split("\n")
    strList[0] += "\n"+topLevelCmt
    return "\n".join(strList)


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

msg = 'generate Xdmf wrapper for HDF5 file.'
parser = OptionParser(option_class=damask.extendableOption,
                      usage='%prog options [file[s]]',
                      description = msg,
                      version = scriptID)

(options, filenames) = parser.parse_args()

h5f = filenames[0]

# ----- parse HDF5 file ----- #
h5f_dataDim = {}
h5f_dataPath = {}
h5f_dataType = {}
with h5py.File(h5f, 'a') as f:
    labels = f.keys()
    labels += f['/Constitutive'].keys()
    labels += f['/Crystallite'].keys()
    labels += ['Vx', 'Vy', "Vz"]
    # remove group names as they do not contain real data
    # TODO: use h5py/H5table API to detect dataset name to
    #       avoid necessary name space pruning.
    labels.remove('Constitutive')
    labels.remove('Crystallite')
    labels.remove('Geometry')
    # loop through remaining labels
    for label in labels:
        dataType, h5Path = lables_to_path(label)
        h5f_dataType[label] = dataType
        h5f_dataDim[label] = " ".join(map(str,f[h5Path].shape))
        h5f_dataPath[label] = h5Paths

# ----- constructing xdmf elements ----- #
root = ET.Element("Xdmf", version='3.3')
root.set('xmlns:xi', "http://www.w3.org/2001/XInclude")
root.append(ET.Comment('Generated Xdmf wapper for DAMASH H5 output'))

# usually there should only be ONE domain
domain = ET.SubElement(root, 'Domain',
                       Name=h5f_base.split(".")[0])

# use global topology through reference
grid = ET.SubElement(domain, 'Grid', GridType="Uniform")
# geometry section
geometry = ET.SubElement(grid, 'Geometry', GeometryType="VXVYVZ")
for vector in ["Vz", "Vy", "Vx"]:
    dataitem = ET.SubElement(geometry, "DataItem",
                             DataType="Float",
                             Dimensions=h5f_dataDim[vector],
                             Name=vector,
                             Format="HDF")
    dataitem.text = h5f_base.split("/")[-1] + ":{}".format(h5f_dataPath[vector])
# topology section
# TODO: support for other format based on given option
meshDim = [h5f_dataDim["Vz"], h5f_dataDim["Vy"], h5f_dataDim["Vx"]]
topology = ET.SubElement(grid, 'Topology',
                         TopologyType="3DRectMesh",
                         Dimensions=" ".join(map(str, meshDim)))

# attributes section
# Question: how to properly handle data mapping for multiphase situations
labelsProcessed = ['Vx', 'Vy', 'Vz']
# walk through each attributes
for label in labels:
    if label in labelsProcessed: continue
    print "adding {}...".format(label)
    attr = ET.SubElement(grid, 'Attribute',
                         Name=label,
                         Type="None",
                         Center="Cell")
    dataitem = ET.SubElement(attr, 'DataItem',
                             Name=label,
                             Format='HDF',
                             Dimensions=h5f_dataDim[label])
    dataitem.text = h5f_base + ":" + h5f_dataPath[label]
    # update progress list
    labelsProcessed.append(label)


# pretty print the xdmf(xml) file content
xmlstr = minidom.parseString(ET.tostring(root)).toprettyxml(indent="\t")
xmlstr = addTopLvlCmt(xmlstr, '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>')
# write str to file through native python API
with open(h5f.replace(".h5", ".xmf"), 'w') as f:
    f.write(xmlstr)
