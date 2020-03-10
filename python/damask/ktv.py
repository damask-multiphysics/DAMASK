import numpy as np
import vtk
#from vtk.util import numpy_support

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
        coordArray = [vtk.vtkDoubleArray(),vtk.vtkDoubleArray(),vtk.vtkDoubleArray()]
        for dim in [0,1,2]:
            for c in np.linspace(0,size[dim],1+grid[dim]):
                coordArray[dim].InsertNextValue(c)

        geom = vtk.vtkRectilinearGrid()
        geom.SetDimensions(*(grid+1))
        geom.SetXCoordinates(coordArray[0])
        geom.SetYCoordinates(coordArray[1])
        geom.SetZCoordinates(coordArray[2])
        
        return VTK(geom)


    @staticmethod
    def from_unstructuredGrid(nodes,connectivity,elem):
        geom = vtk.vtkUnstructuredGrid()
        geom.SetPoints(nodes)
        geom.Allocate(connectivity.shape[0])

        if   elem == 'TRIANGLE':
            vtk_type = vtk.VTK_TRIANGLE
            n_nodes = 3
        elif elem == 'QUAD':
            vtk_type = vtk.VTK_QUAD
            n_nodes = 4
        elif elem == 'TETRA':
            vtk_type = vtk.VTK_TETRA
            n_nodes = 4
        elif elem == 'HEXAHEDRON':
            vtk_type = vtk.VTK_HEXAHEDRON
            n_nodes = 8

        for i in connectivity:
            geom.InsertNextCell(vtk_type,n_nodes,i-1)

        return VTK(geom)

    def write(self,fname):
        print('tbd',fname)
