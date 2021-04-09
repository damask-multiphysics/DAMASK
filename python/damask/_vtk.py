import os
import multiprocessing as mp
from pathlib import Path

import numpy as np
import numpy.ma as ma
import vtk
from vtk.util.numpy_support import numpy_to_vtk            as np_to_vtk
from vtk.util.numpy_support import numpy_to_vtkIdTypeArray as np_to_vtkIdTypeArray
from vtk.util.numpy_support import vtk_to_numpy            as vtk_to_np

from . import util
from . import Table


class VTK:
    """
    Spatial visualization (and potentially manipulation).

    High-level interface to VTK.
    """

    def __init__(self,vtk_data):
        """
        New spatial visualization.

        Parameters
        ----------
        vtk_data : subclass of vtk.vtkDataSet
            Description of geometry and topology, optionally with attached data.
            Valid types are vtk.vtkRectilinearGrid, vtk.vtkUnstructuredGrid,
            or vtk.vtkPolyData.

        """
        self.vtk_data = vtk_data


    @staticmethod
    def from_rectilinear_grid(grid,size,origin=np.zeros(3)):
        """
        Create VTK of type vtk.vtkRectilinearGrid.

        This is the common type for grid solver results.

        Parameters
        ----------
        grid : iterable of int, len (3)
            Number of cells along each dimension.
        size : iterable of float, len (3)
            Physical lengths along each dimension.
        origin : iterable of float, len (3), optional
            Spatial origin coordinates.

        """
        vtk_data = vtk.vtkRectilinearGrid()
        vtk_data.SetDimensions(*(np.array(grid)+1))
        coord = [np_to_vtk(np.linspace(origin[i],origin[i]+size[i],grid[i]+1),deep=True) for i in [0,1,2]]
        [coord[i].SetName(n) for i,n in enumerate(['x','y','z'])]
        vtk_data.SetXCoordinates(coord[0])
        vtk_data.SetYCoordinates(coord[1])
        vtk_data.SetZCoordinates(coord[2])

        return VTK(vtk_data)


    @staticmethod
    def from_unstructured_grid(nodes,connectivity,cell_type):
        """
        Create VTK of type vtk.vtkUnstructuredGrid.

        This is the common type for FEM solver results.

        Parameters
        ----------
        nodes : numpy.ndarray of shape (:,3)
            Spatial position of the nodes.
        connectivity : numpy.ndarray of np.dtype = int
            Cell connectivity (0-based), first dimension determines #Cells,
            second dimension determines #Nodes/Cell.
        cell_type : str
            Name of the vtk.vtkCell subclass. Tested for TRIANGLE, QUAD, TETRA, and HEXAHEDRON.

        """
        vtk_nodes = vtk.vtkPoints()
        vtk_nodes.SetData(np_to_vtk(nodes))
        cells = vtk.vtkCellArray()
        cells.SetNumberOfCells(connectivity.shape[0])
        T = np.concatenate((np.ones((connectivity.shape[0],1),dtype=np.int64)*connectivity.shape[1],
                            connectivity),axis=1).ravel()
        cells.SetCells(connectivity.shape[0],np_to_vtkIdTypeArray(T,deep=True))

        vtk_data = vtk.vtkUnstructuredGrid()
        vtk_data.SetPoints(vtk_nodes)
        cell_types = {'TRIANGLE':vtk.VTK_TRIANGLE, 'QUAD':vtk.VTK_QUAD,
                      'TETRA'   :vtk.VTK_TETRA,    'HEXAHEDRON':vtk.VTK_HEXAHEDRON}
        vtk_data.SetCells(cell_types[cell_type.split("_",1)[-1].upper()],cells)

        return VTK(vtk_data)


    @staticmethod
    def from_poly_data(points):
        """
        Create VTK of type vtk.polyData.

        This is the common type for point-wise data.

        Parameters
        ----------
        points : numpy.ndarray of shape (:,3)
            Spatial position of the points.

        """
        N = points.shape[0]
        vtk_points = vtk.vtkPoints()
        vtk_points.SetData(np_to_vtk(points))

        vtk_cells = vtk.vtkCellArray()
        vtk_cells.SetNumberOfCells(N)
        vtk_cells.SetCells(N,np_to_vtkIdTypeArray(np.stack((np.ones  (N,dtype=np.int64),
                                                            np.arange(N,dtype=np.int64)),axis=1).ravel(),deep=True))

        vtk_data = vtk.vtkPolyData()
        vtk_data.SetPoints(vtk_points)
        vtk_data.SetVerts(vtk_cells)

        return VTK(vtk_data)


    @staticmethod
    def load(fname,dataset_type=None):
        """
        Load from VTK file.

        Parameters
        ----------
        fname : str or pathlib.Path
            Filename for reading. Valid extensions are .vtr, .vtu, .vtp, and .vtk.
        dataset_type : str, optional
            Name of the vtk.vtkDataSet subclass when opening a .vtk file. Valid types are vtkRectilinearGrid,
            vtkUnstructuredGrid, and vtkPolyData.

        """
        if not os.path.isfile(fname):                                                               # vtk has a strange error handling
            raise FileNotFoundError(f'no such file: {fname}')
        ext = Path(fname).suffix
        if ext == '.vtk' or dataset_type is not None:
            reader = vtk.vtkGenericDataObjectReader()
            reader.SetFileName(str(fname))
            if dataset_type is None:
                raise TypeError('Dataset type for *.vtk file not given.')
            elif dataset_type.lower().endswith('rectilineargrid'):
                reader.Update()
                vtk_data = reader.GetRectilinearGridOutput()
            elif dataset_type.lower().endswith('unstructuredgrid'):
                reader.Update()
                vtk_data = reader.GetUnstructuredGridOutput()
            elif dataset_type.lower().endswith('polydata'):
                reader.Update()
                vtk_data = reader.GetPolyDataOutput()
            else:
                raise TypeError(f'Unknown dataset type {dataset_type} for vtk file')
        else:
            if   ext == '.vtr':
                reader = vtk.vtkXMLRectilinearGridReader()
            elif ext == '.vtu':
                reader = vtk.vtkXMLUnstructuredGridReader()
            elif ext == '.vtp':
                reader = vtk.vtkXMLPolyDataReader()
            else:
                raise TypeError(f'Unknown file extension {ext}')

            reader.SetFileName(str(fname))
            reader.Update()
            vtk_data = reader.GetOutput()

        return VTK(vtk_data)


    @staticmethod
    def _write(writer):
        """Wrapper for parallel writing."""
        writer.Write()
    def save(self,fname,parallel=True,compress=True):
        """
        Save as VTK file.

        Parameters
        ----------
        fname : str or pathlib.Path
            Filename for writing.
        parallel : boolean, optional
            Write data in parallel background process. Defaults to True.
        compress : bool, optional
            Compress with zlib algorithm. Defaults to True.

        """
        if   isinstance(self.vtk_data,vtk.vtkRectilinearGrid):
            writer = vtk.vtkXMLRectilinearGridWriter()
        elif isinstance(self.vtk_data,vtk.vtkUnstructuredGrid):
            writer = vtk.vtkXMLUnstructuredGridWriter()
        elif isinstance(self.vtk_data,vtk.vtkPolyData):
            writer = vtk.vtkXMLPolyDataWriter()

        default_ext = '.'+writer.GetDefaultFileExtension()
        ext = Path(fname).suffix
        writer.SetFileName(str(fname)+(default_ext if default_ext != ext else ''))

        if compress:
            writer.SetCompressorTypeToZLib()
        else:
            writer.SetCompressorTypeToNone()
        writer.SetDataModeToBinary()
        writer.SetInputData(self.vtk_data)

        if parallel:
            try:
                mp_writer = mp.Process(target=self._write,args=(writer,))
                mp_writer.start()
            except TypeError:
                writer.Write()
        else:
            writer.Write()


    # Check https://blog.kitware.com/ghost-and-blanking-visibility-changes/ for missing data
    # Needs support for damask.Table
    def add(self,data,label=None):
        """
        Add data to either cells or points.

        Parameters
        ----------
        data : numpy.ndarray or numpy.ma.MaskedArray
            Data to add. First dimension needs to match either
            number of cells or number of points.
        label : str
            Data label.

        """
        N_points = self.vtk_data.GetNumberOfPoints()
        N_cells  = self.vtk_data.GetNumberOfCells()

        if isinstance(data,np.ndarray):
            if label is None:
                raise ValueError('No label defined for numpy.ndarray')

            N_data = data.shape[0]
            data_ = np.where(data.mask,data.fill_value,data) if isinstance(data,ma.MaskedArray) else\
                    data
            d = np_to_vtk((data_.astype(np.single) if data_.dtype in [np.double, np.longdouble] else
                           data_).reshape(N_data,-1),deep=True)                                      # avoid large files
            d.SetName(label)

            if   N_data == N_points:
                self.vtk_data.GetPointData().AddArray(d)
            elif N_data == N_cells:
                self.vtk_data.GetCellData().AddArray(d)
            else:
                raise ValueError(f'Cell / point count ({N_cells} / {N_points}) differs from data ({N_data}).')
        elif isinstance(data,Table):
            raise NotImplementedError('damask.Table')
        else:
            raise TypeError


    def get(self,label):
        """
        Get either cell or point data.

        Cell data takes precedence over point data, i.e. this
        function assumes that labels are unique among cell and
        point data.

        Parameters
        ----------
        label : str
            Data label.

        """
        cell_data = self.vtk_data.GetCellData()
        for a in range(cell_data.GetNumberOfArrays()):
            if cell_data.GetArrayName(a) == label:
                return vtk_to_np(cell_data.GetArray(a))

        point_data = self.vtk_data.GetPointData()
        for a in range(point_data.GetNumberOfArrays()):
            if point_data.GetArrayName(a) == label:
                return vtk_to_np(point_data.GetArray(a))

        raise ValueError(f'Array "{label}" not found.')


    def get_comments(self):
        """Return the comments."""
        fielddata = self.vtk_data.GetFieldData()
        for a in range(fielddata.GetNumberOfArrays()):
            if fielddata.GetArrayName(a) == 'comments':
                comments = fielddata.GetAbstractArray(a)
                return [comments.GetValue(i) for i in range(comments.GetNumberOfValues())]
        return []


    def set_comments(self,comments):
        """
        Set comments.

        Parameters
        ----------
        comments : str or list of str
            Comments.

        """
        s = vtk.vtkStringArray()
        s.SetName('comments')
        for c in [comments] if isinstance(comments,str) else comments:
            s.InsertNextValue(c)
        self.vtk_data.GetFieldData().AddArray(s)


    def add_comments(self,comments):
        """
        Add comments.

        Parameters
        ----------
        comments : str or list of str
            Comments to add.

        """
        self.set_comments(self.get_comments() + ([comments] if isinstance(comments,str) else comments))


    def __repr__(self):
        """ASCII representation of the VTK data."""
        writer = vtk.vtkDataSetWriter()
        writer.SetHeader(f'# {util.execution_stamp("VTK")}')
        writer.WriteToOutputStringOn()
        writer.SetInputData(self.vtk_data)
        writer.Write()
        return writer.GetOutputString()


    def show(self):
        """
        Render.

        See http://compilatrix.com/article/vtk-1 for further ideas.
        """
        try:
            import wx
            _ = wx.App(False)                                                                       # noqa
            width, height = wx.GetDisplaySize()
        except ImportError:
            try:
                import tkinter
                tk = tkinter.Tk()
                width  = tk.winfo_screenwidth()
                height = tk.winfo_screenheight()
                tk.destroy()
            except Exception as e:
                width  = 1024
                height =  768

        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(self.vtk_data)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        ren = vtk.vtkRenderer()

        window = vtk.vtkRenderWindow()
        window.AddRenderer(ren)

        ren.AddActor(actor)
        ren.SetBackground(0.2,0.2,0.2)

        window.SetSize(width,height)

        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(window)

        iren.Initialize()
        window.Render()
        iren.Start()
