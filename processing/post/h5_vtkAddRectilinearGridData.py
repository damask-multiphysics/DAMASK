#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os
import vtk
import damask
from vtk.util import numpy_support
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID = ' '.join([scriptName, damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
msg = "Add scalars, vectors, and/or an RGB tuple from"
msg += "an HDF5 to existing VTK rectilinear grid (.vtr/.vtk)."
parser = OptionParser(option_class=damask.extendableOption,
                      usage='%prog options [file[s]]',
                      description=msg,
                      version=scriptID)
parser.add_option('--vtk',
                  dest='vtk',
                  type='string', metavar='string',
                  help='VTK file name')
parser.add_option('--inplace',
                  dest='inplace',
                  action='store_true',
                  help='modify VTK file in-place')
parser.add_option('-r', '--render',
                  dest='render',
                  action='store_true',
                  help='open output in VTK render window')
parser.add_option('-d', '--data',
                  dest='data',
                  action='extend', metavar='<string LIST>',
                  help='scalar/vector value(s) label(s)')
parser.add_option('-t', '--tensor',
                  dest='tensor',
                  action='extend', metavar='<string LIST>',
                  help='tensor (3x3) value label(s)')
parser.add_option('-c', '--color',
                  dest='color',
                  action='extend', metavar='<string LIST>',
                  help='RGB color tuple label')
parser.add_option('-m',
                  '--mode',
                  dest='mode',
                  metavar='string',
                  type='choice', choices=['cell', 'point'],
                  help='cell-centered or point-centered coordinates')

parser.set_defaults(data=[],
                    tensor=[],
                    color=[],
                    mode='cell',
                    inplace=False,
                    render=False)

(options, filenames) = parser.parse_args()

# ----- Legacy VTK format support ----- #
if os.path.splitext(options.vtk)[1] == '.vtr':
    reader = vtk.vtkXMLRectilinearGridReader()
    reader.SetFileName(options.vtk)
    reader.Update()
    rGrid = reader.GetOutput()
elif os.path.splitext(options.vtk)[1] == '.vtk':
    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(options.vtk)
    reader.Update()
    rGrid = reader.GetRectilinearGridOutput()
else:
    parser.error('Unsupported VTK file type extension.')

Npoints = rGrid.GetNumberOfPoints()
Ncells = rGrid.GetNumberOfCells()

# ----- Summary output (Sanity Check) ----- #
msg = '{}: {} points and {} cells...'.format(options.vtk,
                                             Npoints,
                                             Ncells)
damask.util.croak(msg)

# ----- Read HDF5 file ----- #
# NOTE:
# --> It is possible in the future we are trying to add data
#     from different increment into the same VTK file, but
#     this feature is not supported for the moment.
# --> Let it fail, if the HDF5 is invalid, python interpretor
# --> should be able to catch this error.
h5f = damask.H5Table(filenames[0], new_file=False)

# ----- Process data ----- #
featureToAdd = {'data': options.data,
                'tensor': options.tensor,
                'color': options.color}
VTKarray = {}  # store all vtkData in dict, then ship them to file
for dataType in featureToAdd.keys():
    featureNames = featureToAdd[dataType]
    for featureName in featureNames:
        VTKtype = vtk.VTK_DOUBLE
        VTKdata = h5f.get_data(featureName)
        if dataType == 'color':
            VTKtype = vtk.VTK_UNSIGNED_CHAR
            VTKdata = (VTKdata*255).astype(int)
        elif dataType == 'tensor':
            # Force symmetries tensor type data
            VTKdata[:, 1] = VTKdata[:, 3] = 0.5*(VTKdata[:, 1]+VTKdata[:, 3])
            VTKdata[:, 2] = VTKdata[:, 6] = 0.5*(VTKdata[:, 2]+VTKdata[:, 6])
            VTKdata[:, 5] = VTKdata[:, 7] = 0.5*(VTKdata[:, 5]+VTKdata[:, 7])
        # use vtk build-in numpy support to add data (much faster)
        # NOTE:
        # --> deep copy is necessary here, otherwise memory leak could occur
        VTKarray[featureName] = numpy_support.numpy_to_vtk(num_array=VTKdata,
                                                           deep=True,
                                                           array_type=VTKtype)
        VTKarray[featureName].SetName(featureName)

# ----- ship data to vtkGrid ----- #
mode = options.mode
damask.util.croak('{} mode...'.format(mode))

# NOTE:
# --> For unknown reason, Paraview only recognize one
#     tensor attributes per cell, thus it would be safe
#     to only add one attributes as tensor.
for dataType in featureToAdd.keys():
    featureNames = featureToAdd[dataType]
    for featureName in featureNames:
        if dataType == 'color':
            if mode == 'cell':
                rGrid.GetCellData().SetScalars(VTKarray[featureName])
            elif mode == 'point':
                rGrid.GetPointData().SetScalars(VTKarray[featureName])
        elif dataType == 'tensor':
            if mode == 'cell':
                rGrid.GetCellData().SetTensors(VTKarray[featureName])
            elif mode == 'point':
                rGrid.GetPointData().SetTensors(VTKarray[featureName])
        else:
            if mode == 'cell':
                rGrid.GetCellData().AddArray(VTKarray[featureName])
            elif mode == 'point':
                rGrid.GetPointData().AddArray(VTKarray[featureName])

rGrid.Modified()
if vtk.VTK_MAJOR_VERSION <= 5:
    rGrid.Update()

# ----- write Grid to VTK file ----- #
writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetDataModeToBinary()
writer.SetCompressorTypeToZLib()
vtkFileN = os.path.splitext(options.vtk)[0]
vtkExtsn = '.vtr' if options.inplace else '_added.vtr'
writer.SetFileName(vtkFileN+vtkExtsn)
if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInput(rGrid)
else:
    writer.SetInputData(rGrid)
writer.Write()

# ----- render results from script ----- #
if options.render:
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(rGrid)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Create the graphics structure. The renderer renders into the
    # render window. The render window interactor captures mouse events
    # and will perform appropriate camera or actor manipulation
    # depending on the nature of the events.

    ren = vtk.vtkRenderer()

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)

    ren.AddActor(actor)
    ren.SetBackground(1, 1, 1)
    renWin.SetSize(200, 200)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    iren.Initialize()
    renWin.Render()
    iren.Start()
