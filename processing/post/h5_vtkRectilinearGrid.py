#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

# ------------------------------------------------------------------ #
# NOTE:                                                              #
#   1. It might be a good idea to separate IO and calculation.       #
#   2. Some of the calculation could be useful in other situations,  #
#      why not build a math_util, or math_sup module that contains   #
#      all the useful functions.                                     #
# ------------------------------------------------------------------ #

import os
import vtk
import numpy as np
import damask
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID = ' '.join([scriptName, damask.version])


# ----- HELPER FUNCTION ----- #
def getMeshFromXYZ(xyzArray, mode):
    """calc Vx,Vy,Vz vectors for vtk rectangular mesh"""
    # NOTE:
    # --> np.unique will automatically sort the list
    # --> although not exactly n(1), but since mesh dimension should
    #     small anyway, so this is still light weight.
    dim = xyzArray.shape[1]  # 2D:2, 3D:3
    coords = [np.unique(xyzArray[:, i]) for i in xrange(dim)]

    if mode == 'cell':
        # since x, y, z might now have the same number of elements,
        # we have to deal with them individually
        for ri in xrange(dim):
            vctr_pt = coords[ri]
            vctr_cell = np.empty(len(vctr_pt)+1)
            # calculate first and last end point
            vctr_cell[0] = vctr_pt[0] - 0.5*abs(vctr_pt[1] - vctr_pt[0])
            vctr_cell[-1] = vctr_pt[-1] + 0.5*abs(vctr_pt[-2] - vctr_pt[-1])
            for cj in xrange(1, len(vctr_cell)-1):
                vctr_cell[cj] = 0.5*(vctr_pt[cj-1] + vctr_pt[cj])
            # update the coords
            coords[ri] = vctr_cell

    if dim < 3:
        coords.append([0])  # expand to a 3D with 0 for z

    # auxiliary description
    grid = np.array(map(len, coords), 'i')
    N = grid.prod() if mode == 'point' else (grid-1).prod()
    return coords, grid, N

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

msg = "Create regular voxel grid from points in an ASCIItable."
parser = OptionParser(option_class=damask.extendableOption,
                      usage='%prog options [file[s]]',
                      description=msg,
                      version=scriptID)

parser.add_option('-m',
                  '--mode',
                  dest='mode',
                  metavar='string',
                  type='choice', choices=['cell', 'point'],
                  help='cell-centered or point-centered coordinates')
parser.add_option('-p',
                  '--pos', '--position',
                  dest='pos',
                  type='string', metavar='string',
                  help='label of coordinates [%default]')

parser.set_defaults(mode='cell',
                    pos='pos')

(options, filenames) = parser.parse_args()

# ----- loop over input files ----- #
for name in filenames:
    try:
        h5f = damask.H5Table(name, new_file=False)
    except:
        continue
    damask.util.report(scriptName, name)

    # ----- read xyzArray from HDF5 file ----- #
    xyzArray = h5f.get_data(options.pos)

    # ----- figure out size and grid ----- #
    coords, grid, N = getMeshFromXYZ(xyzArray, options.mode)

    # ----- process data ----- #
    rGrid = vtk.vtkRectilinearGrid()
    # WARNING: list expansion does not work here as these are
    #          just pointers for a vtk instance. Simply put,
    #          DON't USE
    #            [<VTK_CONSTRUCTOR>] * <NUM_OF_ELEMENTS>
    coordArray = [vtk.vtkDoubleArray(),
                  vtk.vtkDoubleArray(),
                  vtk.vtkDoubleArray()]

    rGrid.SetDimensions(*grid)
    for i, points in enumerate(coords):
        for point in points:
            coordArray[i].InsertNextValue(point)

    rGrid.SetXCoordinates(coordArray[0])
    rGrid.SetYCoordinates(coordArray[1])
    rGrid.SetZCoordinates(coordArray[2])

    # ----- output result ----- #
    dirPath = os.path.split(name)[0]
    if name:
        writer = vtk.vtkXMLRectilinearGridWriter()
        writer.SetCompressorTypeToZLib()
        writer.SetDataModeToBinary()
        # getting the name is a little bit tricky
        vtkFileName = os.path.splitext(os.path.split(name)[1])[0]
        vtkFileName += '_{}({})'.format(options.pos, options.mode)
        vtkFileName += '.' + writer.GetDefaultFileExtension()
        writer.SetFileName(os.path.join(dirPath, vtkFileName))
    else:
        writer = vtk.vtkDataSetWriter()
        writer.SetHeader('# powered by '+scriptID)
        writer.WriteToOutputStringOn()

    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(rGrid)
    else:
        writer.SetInputData(rGrid)

    writer.Write()
