import os

import pandas as pd
import numpy as np
import vtk
from vtk.util import numpy_support as nps

from . import table
from . import version

class VTK: # capitals needed/preferred?
    """
    Manage vtk files.

    tbd
    """

    def __init__(self,geom):
        """tbd."""
        self.geom = geom

    @staticmethod
    def from_rectilinearGrid(grid,size,origin=np.zeros(3)):
        """Check https://blog.kitware.com/ghost-and-blanking-visibility-changes/ for missing data."""
        coordArray = [vtk.vtkDoubleArray(),vtk.vtkDoubleArray(),vtk.vtkDoubleArray()]
        for dim in [0,1,2]:
            for c in np.linspace(origin[dim],origin[dim]+size[dim],grid[dim]):
                coordArray[dim].InsertNextValue(c)

        geom = vtk.vtkRectilinearGrid()
        geom.SetDimensions(*grid)
        geom.SetXCoordinates(coordArray[0])
        geom.SetYCoordinates(coordArray[1])
        geom.SetZCoordinates(coordArray[2])

        return VTK(geom)


    @staticmethod
    def from_unstructuredGrid(nodes,connectivity,cell_type):
        """
        Create an unstructured grid (mesh).

        connectivity: 0 based at the moment, shape Ncell x N nodes
        cell_type: TRIANGLE, 'QUAD', 'TETRA','HEXAHEDRON'

        """
        vtk_nodes = vtk.vtkPoints()
        vtk_nodes.SetData(nps.numpy_to_vtk(nodes))
        cells = vtk.vtkCellArray()
        cells.SetNumberOfCells(connectivity.shape[0])
        T = np.concatenate((np.ones((connectivity.shape[0],1),dtype=np.int64)*connectivity.shape[1],
                           connectivity),axis=1).ravel()
        cells.SetCells(connectivity.shape[0],nps.numpy_to_vtk(T, deep=True, array_type=vtk.VTK_ID_TYPE))

        geom = vtk.vtkUnstructuredGrid()
        geom.SetPoints(vtk_nodes)
        geom.SetCells(eval('vtk.VTK_{}'.format(cell_type.upper())),cells)

        return VTK(geom)


    @staticmethod
    def from_polyData(points):
        vtk_points= vtk.vtkPoints()
        vtk_points.SetData(nps.numpy_to_vtk(points))

        vertices = vtk.vtkCellArray()
        vertices.SetNumberOfCells(points.shape[0])
        T = np.concatenate((np.ones((points.shape[0],1),dtype=np.int64),
                           np.arange(points.shape[0],dtype=np.int64).reshape(-1,1)),axis=1).ravel()
        vertices.SetCells(points.shape[0],nps.numpy_to_vtk(T, deep=True, array_type=vtk.VTK_ID_TYPE))

        geom = vtk.vtkPolyData()
        geom.SetPoints(vtk_points)
        geom.SetVerts(vertices)

        return VTK(geom)

    @staticmethod
    def from_file(fname,ftype=None):
        ext = os.path.splitext(fname)[1]
        if ext == '.vtk':
            reader = vtk.vtkGenericDataObjectReader()
            reader.SetFileName(fname)
            reader.Update()
            if   ftype.lower() == 'rectilineargrid':
                geom = reader.GetRectilinearGridOutput()
            elif ftype.lower() == 'unstructuredgrid':
                geom = reader.GetUnstructuredGridOutput()
            elif ftype.lower() == 'polydata':
                geom = reader.GetPolyDataOutput()
            else:
                raise Exception
        else:
            if   ext == '.vtr':
                reader = vtk.vtkXMLRectilinearGridReader()
            elif ext == '.vtu':
                reader = vtk.vtkXMLUnstructuredGridReader()
            elif ext == '.vtp':
                reader = vtk.vtkXMLPolyDataReader()
            else:
                raise Exception

            reader.SetFileName(fname)
            reader.Update()
            geom = reader.GetOutput()

        return VTK(geom)


    def write(self,fname):
        """ToDo: Check if given fileextension makes sense."""
        if  (isinstance(self.geom,vtk.vtkRectilinearGrid)):
            writer = vtk.vtkXMLRectilinearGridWriter()
        elif(isinstance(self.geom,vtk.vtkUnstructuredGrid)):
            writer = vtk.vtkXMLUnstructuredGridWriter()
        elif(isinstance(self.geom,vtk.vtkPolyData)):
            writer = vtk.vtkXMLPolyDataWriter()

        writer.SetFileName('{}.{}'.format(os.path.splitext(fname)[0],
                                          writer.GetDefaultFileExtension()))
        writer.SetCompressorTypeToZLib()
        writer.SetDataModeToBinary()
        writer.SetInputData(self.geom)

        writer.Write()


    def add(self,data,label=None):

        Npoints = self.geom.GetNumberOfPoints()
        Ncells  = self.geom.GetNumberOfCells()

        if   isinstance(data,np.ndarray):
            shape = [data.shape[0],np.product(data.shape[1:],dtype=np.int)]
            d = nps.numpy_to_vtk(num_array=data.reshape(shape),deep=True)
            d.SetName(label)
            if   shape[0] == Ncells:
                self.geom.GetCellData().AddArray(d)
            elif shape[0] == Npoints:
                self.geom.GetPointData().AddArray(d)
        elif isinstance(data,pd.DataFrame):
            pass
        elif isinstance(data,table):
            pass


    def __repr__(self):
        """ASCII representation of the VTK data."""
        writer = vtk.vtkDataSetWriter()
        writer.SetHeader('# DAMASK.VTK v{}'.format(version))
        writer.WriteToOutputStringOn()
        writer.SetInputData(self.geom)
        writer.Write()
        return writer.GetOutputString()
