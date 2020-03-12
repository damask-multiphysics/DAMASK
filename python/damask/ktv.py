import os

import pandas as pd
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk as np_to_vtk

from . import table
from . import version

class VTK:
    """
    Spatial visualization (and potentially manipulation).

    High-level interface to VTK.
    """

    def __init__(self,geom):
        """
        Set geometry and topology.

        Parameters
        ----------
        geom : subclass of vtk.vtkDataSet
          Description of geometry and topology. Valid types are vtk.vtkRectilinearGrid,
          vtk.vtkUnstructuredGrid, or vtk.vtkPolyData.

        """
        self.geom = geom

    @staticmethod
    def from_rectilinearGrid(grid,size,origin=np.zeros(3)):
        """
        Create VTK of type vtk.vtkRectilinearGrid.

        This is the common type for results from the grid solver.

        Parameters
        ----------
        grid : numpy.ndarray of shape (3) of np.dtype = int
          Number of cells.
        size : numpy.ndarray of shape (3)
          Physical length.
        origin : numpy.ndarray of shape (3), optional
          Spatial origin.

        """
        coordArray = [vtk.vtkDoubleArray(),vtk.vtkDoubleArray(),vtk.vtkDoubleArray()]
        for dim in [0,1,2]:
            coords = np.linspace(origin[dim],origin[dim]+size[dim],grid[dim]+1)
            coordArray[dim].SetArray(np_to_vtk(coords),grid[dim]+1,1)

        geom = vtk.vtkRectilinearGrid()
        geom.SetDimensions(*(grid+1))
        geom.SetXCoordinates(coordArray[0])
        geom.SetYCoordinates(coordArray[1])
        geom.SetZCoordinates(coordArray[2])

        return VTK(geom)


    @staticmethod
    def from_unstructuredGrid(nodes,connectivity,cell_type):
        """
        Create VTK of type vtk.vtkUnstructuredGrid.

        This is the common type for results from FEM solvers.

        Parameters
        ----------
        nodes : numpy.ndarray of shape (:,3)
          Spatial position of the nodes.
        connectivity : numpy.ndarray of np.dtype = int
          Cell connectivity (0-based), first dimension determines #Cells, second dimension determines #Nodes/Cell.
        cell_type : str
          Name of the vtk.vtkCell subclass. Tested for TRIANGLE, QUAD, and HEXAHEDRON.

        """
        vtk_nodes = vtk.vtkPoints()
        vtk_nodes.SetData(np_to_vtk(nodes))
        cells = vtk.vtkCellArray()
        cells.SetNumberOfCells(connectivity.shape[0])
        T = np.concatenate((np.ones((connectivity.shape[0],1),dtype=np.int64)*connectivity.shape[1],
                            connectivity),axis=1).ravel()
        cells.SetCells(connectivity.shape[0],np_to_vtk(T, deep=True, array_type=vtk.VTK_ID_TYPE))

        geom = vtk.vtkUnstructuredGrid()
        geom.SetPoints(vtk_nodes)
        geom.SetCells(eval('vtk.VTK_{}'.format(cell_type.split('_',1)[-1].upper())),cells)

        return VTK(geom)


    @staticmethod
    def from_polyData(points):
        """
        Create VTK of type vtk.polyData.

        This is the common type for point-wise data.

        Parameters
        ----------
        points : numpy.ndarray of shape (:,3)
          Spatial position of the points.

        """
        vtk_points= vtk.vtkPoints()
        vtk_points.SetData(np_to_vtk(points))

        geom = vtk.vtkPolyData()
        geom.SetPoints(vtk_points)

        return VTK(geom)


    @staticmethod
    def from_file(fname,dataset_type=None):
        """
        Create VTK from file.

        Parameters
        ----------
        fname : str
          Filename for reading. Valid extensions are .vtk, .vtr, .vtu, and .vtp.
        dataset_type : str, optional
          Name of the vtk.vtkDataSet subclass when opening an .vtk file. Valid types are vtkRectilinearGrid,
          vtkUnstructuredGrid, and vtkPolyData.

        """
        ext = os.path.splitext(fname)[1]
        if ext == '.vtk':
            reader = vtk.vtkGenericDataObjectReader()
            reader.SetFileName(fname)
            reader.Update()
            if 'rectilineargrid' in    dataset_type.lower():
                geom = reader.GetRectilinearGridOutput()
            elif 'unstructuredgrid' in dataset_type.lower():
                geom = reader.GetUnstructuredGridOutput()
            elif 'polydata' in         dataset_type.lower():
                geom = reader.GetPolyDataOutput()
            else:
                raise TypeError('Unknown dataset type for vtk file {}'.format(dataset_type))
        else:
            if   ext == '.vtr':
                reader = vtk.vtkXMLRectilinearGridReader()
            elif ext == '.vtu':
                reader = vtk.vtkXMLUnstructuredGridReader()
            elif ext == '.vtp':
                reader = vtk.vtkXMLPolyDataReader()
            else:
                raise TypeError('Unknown file extension {}'.format(ext))

            reader.SetFileName(fname)
            reader.Update()
            geom = reader.GetOutput()

        return VTK(geom)


    # ToDo: If extension is given, check for consistency.
    def write(self,fname):
        """
        Write to file.

        Parameters
        ----------
        fname : str
          Filename for writing.

        """
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


    # Check https://blog.kitware.com/ghost-and-blanking-visibility-changes/ for missing data
    # Needs support for pd.DataFrame and/or table
    def add(self,data,label=None):
        """Add data to either cells or points."""
        N_points = self.geom.GetNumberOfPoints()
        N_cells  = self.geom.GetNumberOfCells()

        if   isinstance(data,np.ndarray):
            d = np_to_vtk(num_array=data.reshape(data.shape[0],-1),deep=True)
            d.SetName(label)
            if   data.shape[0] == N_cells:
                self.geom.GetCellData().AddArray(d)
            elif data.shape[0] == N_points:
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
