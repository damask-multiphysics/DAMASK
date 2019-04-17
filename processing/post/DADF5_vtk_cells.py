#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os,vtk
import numpy as np
import argparse
import damask
from vtk.util import numpy_support

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
parser = argparse.ArgumentParser()

#ToDo:  We need to decide on a way of handling arguments of variable lentght
#https://stackoverflow.com/questions/15459997/passing-integer-lists-to-python

#parser.add_argument('--version', action='version', version='%(prog)s {}'.format(scriptID))
parser.add_argument('filenames', nargs='+',
                    help='DADF5 files')

options = parser.parse_args()


# --- loop over input files ------------------------------------------------------------------------

for filename in options.filenames:
  data = damask.DADF5(filename)
  
  if data.structured:                                                                               # for grid solvers use rectilinear grid
    rGrid = vtk.vtkRectilinearGrid()
    coordArray = [vtk.vtkDoubleArray(),
                  vtk.vtkDoubleArray(),
                  vtk.vtkDoubleArray(),
                 ]

    rGrid.SetDimensions(*(data.grid+1))
    for dim in [0,1,2]:
      for c in np.linspace(0,data.size[dim],1+data.grid[dim]):
        coordArray[dim].InsertNextValue(c)

    rGrid.SetXCoordinates(coordArray[0])
    rGrid.SetYCoordinates(coordArray[1])
    rGrid.SetZCoordinates(coordArray[2])


  for i,inc in enumerate(data.increments):
    data.active['increments'] = [inc]
    x = data.get_dataset_location('xi_sl')[0]
    VTKarray = numpy_support.numpy_to_vtk(num_array=data.read_dataset(x,0),deep=True,array_type= vtk.VTK_DOUBLE)
    VTKarray.SetName('xi_sl')
    rGrid.GetCellData().AddArray(VTKarray)
    if data.structured:
      writer = vtk.vtkXMLRectilinearGridWriter()

    writer.SetCompressorTypeToZLib()
    writer.SetDataModeToBinary()
    writer.SetFileName(os.path.join(os.path.split(filename)[0],
                                    os.path.splitext(os.path.split(filename)[1])[0] +
                                    '_inc{:04d}'.format(i) +                                        # ToDo: adjust to lenght of increments
                                    '.' + writer.GetDefaultFileExtension()))
    if data.structured:
      writer.SetInputData(rGrid)
    
    writer.Write()
