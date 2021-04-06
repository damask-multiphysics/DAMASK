import copy
import multiprocessing as mp
from functools import partial
import os
import warnings

import numpy as np
import pandas as pd
import h5py
from scipy import ndimage, spatial
from vtk.util.numpy_support import vtk_to_numpy as vtk_to_np

from . import VTK
from . import util
from . import grid_filters
from . import Rotation


class Grid:
    """
    Geometry definition for grid solvers.

    Create and manipulate geometry definitions for storage as VTK
    rectiliear grid files ('.vtr' extension). A grid contains the
    material ID (referring to the entry in 'material.yaml') and
    the physical size.
    """

    def __init__(self,material,size,origin=[0.0,0.0,0.0],comments=[]):
        """
        New geometry definition for grid solvers.

        Parameters
        ----------
        material : numpy.ndarray
            Material index array (3D).
        size : list or numpy.ndarray
            Physical size of the grid in meter.
        origin : list or numpy.ndarray, optional
            Physical origin of the grid in meter.
        comments : list of str, optional
            Comment lines.

        """
        self.material = material
        self.size     = size
        self.origin   = origin
        self.comments = comments


    def __repr__(self):
        """Basic information on grid definition."""
        mat_min = np.nanmin(self.material)
        mat_max = np.nanmax(self.material)
        mat_N   = self.N_materials
        return util.srepr([
               f'cells  a b c: {util.srepr(self.cells, " x ")}',
               f'size   x y z: {util.srepr(self.size,  " x ")}',
               f'origin x y z: {util.srepr(self.origin,"   ")}',
               f'# materials: {mat_N}' + ('' if mat_min == 0 and mat_max+1 == mat_N else
                                          f' (min: {mat_min}, max: {mat_max})')
              ])


    def __copy__(self):
        """Create deep copy."""
        return copy.deepcopy(self)

    copy = __copy__


    def __eq__(self,other):
        """
        Test equality of other.

        Parameters
        ----------
        other : damask.Grid
            Grid to compare self against.

        """
        return (np.allclose(other.size,self.size)
            and np.allclose(other.origin,self.origin)
            and np.all(other.cells == self.cells)
            and np.all(other.material == self.material))


    @property
    def material(self):
        """Material indices."""
        return self._material

    @material.setter
    def material(self,material):
        if len(material.shape) != 3:
            raise ValueError(f'invalid material shape {material.shape}')
        elif material.dtype not in np.sctypes['float'] + np.sctypes['int']:
            raise TypeError(f'invalid material data type {material.dtype}')
        else:
            self._material = np.copy(material)

            if self.material.dtype in np.sctypes['float'] and \
               np.all(self.material == self.material.astype(int).astype(float)):
                self._material = self.material.astype(int)


    @property
    def size(self):
        """Physical size of grid in meter."""
        return self._size

    @size.setter
    def size(self,size):
        if len(size) != 3 or any(np.array(size) < 0):
            raise ValueError(f'invalid size {size}')
        else:
            self._size = np.array(size)

    @property
    def origin(self):
        """Coordinates of grid origin in meter."""
        return self._origin

    @origin.setter
    def origin(self,origin):
        if len(origin) != 3:
            raise ValueError(f'invalid origin {origin}')
        else:
            self._origin = np.array(origin)

    @property
    def comments(self):
        """Comments, e.g. history of operations."""
        return self._comments

    @comments.setter
    def comments(self,comments):
        self._comments = [str(c) for c in comments] if isinstance(comments,list) else [str(comments)]


    @property
    def cells(self):
        """Number of cells in x,y,z direction."""
        return np.asarray(self.material.shape)


    @property
    def N_materials(self):
        """Number of (unique) material indices within grid."""
        return np.unique(self.material).size


    @staticmethod
    def load(fname):
        """
        Load from VTK rectilinear grid file.

        Parameters
        ----------
        fname : str or or pathlib.Path
            Grid file to read. Valid extension is .vtr, which will be appended
            if not given.

        """
        v = VTK.load(fname if str(fname).endswith('.vtr') else str(fname)+'.vtr')
        comments = v.get_comments()
        cells = np.array(v.vtk_data.GetDimensions())-1
        bbox  = np.array(v.vtk_data.GetBounds()).reshape(3,2).T

        for i,c in enumerate([v.vtk_data.GetXCoordinates(),v.vtk_data.GetYCoordinates(),v.vtk_data.GetZCoordinates()]):
            if not np.allclose(vtk_to_np(c),np.linspace(bbox[0][i],bbox[1][i],cells[i]+1)):
                raise ValueError('regular grid spacing violated')

        return Grid(material = v.get('material').reshape(cells,order='F'),
                    size = bbox[1] - bbox[0],
                    origin = bbox[0],
                    comments=comments)


    @staticmethod
    def load_ASCII(fname):
        """
        Load from geom file.

        Storing geometry files in ASCII format is deprecated.
        This function will be removed in a future version of DAMASK.

        Parameters
        ----------
        fname : str, pathlib.Path, or file handle
            Geometry file to read.

        """
        warnings.warn('Support for ASCII-based geom format will be removed in DAMASK 3.1.0', DeprecationWarning,2)
        try:
            f = open(fname)
        except TypeError:
            f = fname

        f.seek(0)
        try:
            header_length,keyword = f.readline().split()[:2]
            header_length = int(header_length)
        except ValueError:
            header_length,keyword = (-1, 'invalid')
        if not keyword.startswith('head') or header_length < 3:
            raise TypeError('header length information missing or invalid')

        comments = []
        content = f.readlines()
        for i,line in enumerate(content[:header_length]):
            items = line.split('#')[0].lower().strip().split()
            key = items[0] if items else ''
            if   key == 'grid':
                cells  = np.array([  int(dict(zip(items[1::2],items[2::2]))[i]) for i in ['a','b','c']])
            elif key == 'size':
                size   = np.array([float(dict(zip(items[1::2],items[2::2]))[i]) for i in ['x','y','z']])
            elif key == 'origin':
                origin = np.array([float(dict(zip(items[1::2],items[2::2]))[i]) for i in ['x','y','z']])
            else:
                comments.append(line.strip())

        material = np.empty(cells.prod())                                                           # initialize as flat array
        i = 0
        for line in content[header_length:]:
            items = line.split('#')[0].split()
            if len(items) == 3:
                if   items[1].lower() == 'of':
                    items = np.ones(int(items[0]))*float(items[2])
                elif items[1].lower() == 'to':
                    items = np.linspace(int(items[0]),int(items[2]),
                                        abs(int(items[2])-int(items[0]))+1,dtype=float)
                else:                        items = list(map(float,items))
            else:                            items = list(map(float,items))
            material[i:i+len(items)] = items
            i += len(items)

        if i != cells.prod():
            raise TypeError(f'invalid file: expected {cells.prod()} entries, found {i}')

        if not np.any(np.mod(material,1) != 0.0):                                                   # no float present
            material = material.astype('int') - (1 if material.min() > 0 else 0)

        return Grid(material.reshape(cells,order='F'),size,origin,comments)


    @staticmethod
    def load_DREAM3D(fname,
                     feature_IDs=None,cell_data=None,
                     phases='Phases',Euler_angles='EulerAngles',
                     base_group=None):
        """
        Load DREAM.3D (HDF5) file.

        Data in DREAM.3D files can be stored per cell ('CellData') and/or
        per grain ('Grain Data'). Per default, cell-wise data is assumed.

        damask.ConfigMaterial.load_DREAM3D gives the corresponding material definition.

        Parameters
        ----------
        fname : str
            Filename of the DREAM.3D (HDF5) file.
        feature_IDs : str
            Name of the dataset containing the mapping between cells and
            grain-wise data. Defaults to 'None', in which case cell-wise
            data is used.
        cell_data : str
            Name of the group (folder) containing cell-wise data. Defaults to
            None in wich case it is automatically detected.
        phases : str
            Name of the dataset containing the phase ID. It is not used for
            grain-wise data, i.e. when feature_IDs is not None.
            Defaults to 'Phases'.
        Euler_angles : str
            Name of the dataset containing the crystallographic orientation as
            Euler angles in radians It is not used for grain-wise data, i.e.
            when feature_IDs is not None. Defaults to 'EulerAngles'.
        base_group : str
            Path to the group (folder) that contains geometry (_SIMPL_GEOMETRY),
            and grain- or cell-wise data. Defaults to None, in which case
            it is set as the path that contains _SIMPL_GEOMETRY/SPACING.


        """
        b = util.DREAM3D_base_group(fname)      if base_group is None else base_group
        c = util.DREAM3D_cell_data_group(fname) if cell_data  is None else cell_data
        f = h5py.File(fname, 'r')

        cells  = f['/'.join([b,'_SIMPL_GEOMETRY','DIMENSIONS'])][()]
        size   = f['/'.join([b,'_SIMPL_GEOMETRY','SPACING'])] * cells
        origin = f['/'.join([b,'_SIMPL_GEOMETRY','ORIGIN'])][()]

        if feature_IDs is None:
            phase = f['/'.join([b,c,phases])][()].reshape(-1,1)
            O = Rotation.from_Euler_angles(f['/'.join([b,c,Euler_angles])]).as_quaternion().reshape(-1,4) # noqa
            unique,unique_inverse = np.unique(np.hstack([O,phase]),return_inverse=True,axis=0)
            ma = np.arange(cells.prod()) if len(unique) == cells.prod() else \
                 np.arange(unique.size)[np.argsort(pd.unique(unique_inverse))][unique_inverse]
        else:
            ma = f['/'.join([b,c,feature_IDs])][()].flatten()

        return Grid(ma.reshape(cells,order='F'),size,origin,util.execution_stamp('Grid','load_DREAM3D'))


    @staticmethod
    def from_table(table,coordinates,labels):
        """
        Generate grid from ASCII table.

        Parameters
        ----------
        table : damask.Table
            Table that contains material information.
        coordinates : str
            Label of the vector column containing the spatial coordinates.
            Need to be ordered (1./x fast, 3./z slow).
        labels : str or list of str
            Label(s) of the columns containing the material definition.
            Each unique combination of values results in one material ID.

        """
        cells,size,origin = grid_filters.cellsSizeOrigin_coordinates0_point(table.get(coordinates))

        labels_ = [labels] if isinstance(labels,str) else labels
        unique,unique_inverse = np.unique(np.hstack([table.get(l) for l in labels_]),return_inverse=True,axis=0)

        ma = np.arange(cells.prod()) if len(unique) == cells.prod() else \
             np.arange(unique.size)[np.argsort(pd.unique(unique_inverse))][unique_inverse]

        return Grid(ma.reshape(cells,order='F'),size,origin,util.execution_stamp('Grid','from_table'))


    @staticmethod
    def _find_closest_seed(seeds, weights, point):
        return np.argmin(np.sum((np.broadcast_to(point,(len(seeds),3))-seeds)**2,axis=1) - weights)

    @staticmethod
    def from_Laguerre_tessellation(cells,size,seeds,weights,material=None,periodic=True):
        """
        Generate grid from Laguerre tessellation.

        Parameters
        ----------
        cells : int numpy.ndarray of shape (3)
            Number of cells in x,y,z direction.
        size : list or numpy.ndarray of shape (3)
            Physical size of the grid in meter.
        seeds : numpy.ndarray of shape (:,3)
            Position of the seed points in meter. All points need to lay within the box.
        weights : numpy.ndarray of shape (seeds.shape[0])
            Weights of the seeds. Setting all weights to 1.0 gives a standard Voronoi tessellation.
        material : numpy.ndarray of shape (seeds.shape[0]), optional
            Material ID of the seeds.
            Defaults to None, in which case materials are consecutively numbered.
        periodic : Boolean, optional
            Assume grid to be periodic. Defaults to True.

        """
        if periodic:
            weights_p = np.tile(weights,27)                                                         # Laguerre weights (1,2,3,1,2,3,...,1,2,3)
            seeds_p = np.vstack((seeds  -np.array([size[0],0.,0.]),seeds,  seeds  +np.array([size[0],0.,0.])))
            seeds_p = np.vstack((seeds_p-np.array([0.,size[1],0.]),seeds_p,seeds_p+np.array([0.,size[1],0.])))
            seeds_p = np.vstack((seeds_p-np.array([0.,0.,size[2]]),seeds_p,seeds_p+np.array([0.,0.,size[2]])))
            coords  = grid_filters.coordinates0_point(cells*3,size*3,-size).reshape(-1,3)
        else:
            weights_p = weights
            seeds_p   = seeds
            coords    = grid_filters.coordinates0_point(cells,size).reshape(-1,3)

        pool = mp.Pool(int(os.environ.get('OMP_NUM_THREADS',1)))
        result = pool.map_async(partial(Grid._find_closest_seed,seeds_p,weights_p), [coord for coord in coords])
        pool.close()
        pool.join()
        material_ = np.array(result.get())

        if periodic:
            material_ = material_.reshape(cells*3)
            material_ = material_[cells[0]:cells[0]*2,cells[1]:cells[1]*2,cells[2]:cells[2]*2]%seeds.shape[0]
        else:
            material_ = material_.reshape(cells)

        return Grid(material = material_ if material is None else material[material_],
                    size     = size,
                    comments = util.execution_stamp('Grid','from_Laguerre_tessellation'),
                   )


    @staticmethod
    def from_Voronoi_tessellation(cells,size,seeds,material=None,periodic=True):
        """
        Generate grid from Voronoi tessellation.

        Parameters
        ----------
        cells : int numpy.ndarray of shape (3)
            Number of cells in x,y,z direction.
        size : list or numpy.ndarray of shape (3)
            Physical size of the grid in meter.
        seeds : numpy.ndarray of shape (:,3)
            Position of the seed points in meter. All points need to lay within the box.
        material : numpy.ndarray of shape (seeds.shape[0]), optional
            Material ID of the seeds.
            Defaults to None, in which case materials are consecutively numbered.
        periodic : Boolean, optional
            Assume grid to be periodic. Defaults to True.

        """
        coords = grid_filters.coordinates0_point(cells,size).reshape(-1,3)
        KDTree = spatial.cKDTree(seeds,boxsize=size) if periodic else spatial.cKDTree(seeds)
        devNull,material_ = KDTree.query(coords)

        return Grid(material = (material_ if material is None else material[material_]).reshape(cells),
                    size     = size,
                    comments = util.execution_stamp('Grid','from_Voronoi_tessellation'),
                   )


    _minimal_surface = \
        {'Schwarz P':        lambda x,y,z: np.cos(x) + np.cos(y) + np.cos(z),
         'Double Primitive': lambda x,y,z: ( 0.5 * (np.cos(x)*np.cos(y) + np.cos(y)*np.cos(z) + np.cos(z)*np.cos(x))
                                           + 0.2 * (np.cos(2*x) + np.cos(2*y) + np.cos(2*z)) ),
         'Schwarz D':        lambda x,y,z: ( np.sin(x)*np.sin(y)*np.sin(z)
                                           + np.sin(x)*np.cos(y)*np.cos(z)
                                           + np.cos(x)*np.cos(y)*np.sin(z)
                                           + np.cos(x)*np.sin(y)*np.cos(z) ),
         'Complementary D':  lambda x,y,z: ( np.cos(3*x+y)*np.cos(z) - np.sin(3*x-y)*np.sin(z) + np.cos(x+3*y)*np.cos(z)
                                           + np.sin(x-3*y)*np.sin(z) + np.cos(x-y)*np.cos(3*z) - np.sin(x+y)*np.sin(3*z) ),
         'Double Diamond':   lambda x,y,z: 0.5 * (np.sin(x)*np.sin(y)
                                                + np.sin(y)*np.sin(z)
                                                + np.sin(z)*np.sin(x)
                                                + np.cos(x) * np.cos(y) * np.cos(z) ),
         'Dprime':           lambda x,y,z: 0.5 * ( np.cos(x)*np.cos(y)*np.cos(z)
                                                 + np.cos(x)*np.sin(y)*np.sin(z)
                                                 + np.sin(x)*np.cos(y)*np.sin(z)
                                                 + np.sin(x)*np.sin(y)*np.cos(z)
                                                 - np.sin(2*x)*np.sin(2*y)
                                                 - np.sin(2*y)*np.sin(2*z)
                                                 - np.sin(2*z)*np.sin(2*x) ) - 0.2,
         'Gyroid':           lambda x,y,z: np.cos(x)*np.sin(y) + np.cos(y)*np.sin(z) + np.cos(z)*np.sin(x),
         'Gprime':           lambda x,y,z : ( np.sin(2*x)*np.cos(y)*np.sin(z)
                                            + np.sin(2*y)*np.cos(z)*np.sin(x)
                                            + np.sin(2*z)*np.cos(x)*np.sin(y) ) + 0.32,
         'Karcher K':        lambda x,y,z: ( 0.3 * ( np.cos(x)           + np.cos(y)           + np.cos(z)
                                                   + np.cos(x)*np.cos(y) + np.cos(y)*np.cos(z) + np.cos(z)*np.cos(x) )
                                           - 0.4 * ( np.cos(2*x)         + np.cos(2*y)         + np.cos(2*z) ) ) + 0.2,
         'Lidinoid':         lambda x,y,z: 0.5 * ( np.sin(2*x)*np.cos(y)*np.sin(z)
                                                 + np.sin(2*y)*np.cos(z)*np.sin(x)
                                                 + np.sin(2*z)*np.cos(x)*np.sin(y)
                                                 - np.cos(2*x)*np.cos(2*y)
                                                 - np.cos(2*y)*np.cos(2*z)
                                                 - np.cos(2*z)*np.cos(2*x) ) + 0.15,
         'Neovius':          lambda x,y,z: ( 3 * (np.cos(x)+np.cos(y)+np.cos(z))
                                           + 4 * np.cos(x)*np.cos(y)*np.cos(z) ),
         'Fisher-Koch S':    lambda x,y,z: ( np.cos(2*x)*np.sin(  y)*np.cos(  z)
                                           + np.cos(  x)*np.cos(2*y)*np.sin(  z)
                                           + np.sin(  x)*np.cos(  y)*np.cos(2*z) ),
      }


    @staticmethod
    def from_minimal_surface(cells,size,surface,threshold=0.0,periods=1,materials=(0,1)):
        """
        Generate grid from definition of triply periodic minimal surface.

        Parameters
        ----------
        cells : int numpy.ndarray of shape (3)
            Number of cells in x,y,z direction.
        size : list or numpy.ndarray of shape (3)
            Physical size of the grid in meter.
        surface : str
            Type of the minimal surface. See notes for details.
        threshold : float, optional.
            Threshold of the minimal surface. Defaults to 0.0.
        periods : integer, optional.
            Number of periods per unit cell. Defaults to 1.
        materials : (int, int), optional
            Material IDs. Defaults to (1,2).

        Notes
        -----
        The following triply-periodic minimal surfaces are implemented:
          - Schwarz P
          - Double Primitive
          - Schwarz D
          - Complementary D
          - Double Diamond
          - Dprime
          - Gyroid
          - Gprime
          - Karcher K
          - Lidinoid
          - Neovius
          - Fisher-Koch S

        References
        ----------
        S.B.G. Blanquer et al., Biofabrication 9(2):025001, 2017
        https://doi.org/10.1088/1758-5090/aa6553

        M. Wohlgemuth et al., Macromolecules 34(17):6083-6089, 2001
        https://doi.org/10.1021/ma0019499

        M.-T. Hsieh and L. Valdevit, Software Impacts 6:100026, 2020
        https://doi.org/10.1016/j.simpa.2020.100026

        """
        x,y,z = np.meshgrid(periods*2.0*np.pi*(np.arange(cells[0])+0.5)/cells[0],
                            periods*2.0*np.pi*(np.arange(cells[1])+0.5)/cells[1],
                            periods*2.0*np.pi*(np.arange(cells[2])+0.5)/cells[2],
                            indexing='ij',sparse=True)
        return Grid(material = np.where(threshold < Grid._minimal_surface[surface](x,y,z),materials[1],materials[0]),
                    size     = size,
                    comments = util.execution_stamp('Grid','from_minimal_surface'),
                   )


    def save(self,fname,compress=True):
        """
        Save as VTK rectilinear grid file.

        Parameters
        ----------
        fname : str or pathlib.Path
            Filename to write. Valid extension is .vtr, it will be appended if not given.
        compress : bool, optional
            Compress with zlib algorithm. Defaults to True.

        """
        v = VTK.from_rectilinear_grid(self.cells,self.size,self.origin)
        v.add(self.material.flatten(order='F'),'material')
        v.add_comments(self.comments)

        v.save(fname if str(fname).endswith('.vtr') else str(fname)+'.vtr',parallel=False,compress=compress)


    def save_ASCII(self,fname):
        """
        Save as geom file.

        Storing geometry files in ASCII format is deprecated.
        This function will be removed in a future version of DAMASK.

        Parameters
        ----------
        fname : str or file handle
            Geometry file to write with extension '.geom'.
        compress : bool, optional
            Compress geometry with 'x of y' and 'a to b'.

        """
        warnings.warn('Support for ASCII-based geom format will be removed in DAMASK 3.1.0', DeprecationWarning,2)
        header =  [f'{len(self.comments)+4} header'] + self.comments \
                + ['grid   a {} b {} c {}'.format(*self.cells),
                   'size   x {} y {} z {}'.format(*self.size),
                   'origin x {} y {} z {}'.format(*self.origin),
                   'homogenization 1',
                  ]

        format_string = '%g' if self.material.dtype in np.sctypes['float'] else \
                        '%{}i'.format(1+int(np.floor(np.log10(np.nanmax(self.material)))))
        np.savetxt(fname,
                   self.material.reshape([self.cells[0],np.prod(self.cells[1:])],order='F').T,
                   header='\n'.join(header), fmt=format_string, comments='')


    def show(self):
        """Show on screen."""
        VTK.from_rectilinear_grid(self.cells,self.size,self.origin).show()


    def add_primitive(self,dimension,center,exponent,
                      fill=None,R=Rotation(),inverse=False,periodic=True):
        """
        Insert a primitive geometric object at a given position.

        Parameters
        ----------
        dimension : int or float numpy.ndarray of shape (3)
            Dimension (diameter/side length) of the primitive. If given as
            integers, cell centers are addressed.
            If given as floats, coordinates are addressed.
        center : int or float numpy.ndarray of shape (3)
            Center of the primitive. If given as integers, cell centers are addressed.
            If given as floats, coordinates in space are addressed.
        exponent : numpy.ndarray of shape (3) or float
            Exponents for the three axes.
            0 gives octahedron (ǀxǀ^(2^0) + ǀyǀ^(2^0) + ǀzǀ^(2^0) < 1)
            1 gives sphere     (ǀxǀ^(2^1) + ǀyǀ^(2^1) + ǀzǀ^(2^1) < 1)
        fill : int, optional
            Fill value for primitive. Defaults to material.max()+1.
        R : damask.Rotation, optional
            Rotation of primitive. Defaults to no rotation.
        inverse : Boolean, optional
            Retain original materials within primitive and fill outside.
            Defaults to False.
        periodic : Boolean, optional
            Assume grid to be periodic. Defaults to True.

        """
        # radius and center
        r = np.array(dimension)/2.0*self.size/self.cells if np.array(dimension).dtype in np.sctypes['int'] else \
            np.array(dimension)/2.0
        c = (np.array(center) + .5)*self.size/self.cells if np.array(center).dtype    in np.sctypes['int'] else \
            (np.array(center) - self.origin)

        coords = grid_filters.coordinates0_point(self.cells,self.size,
                                          -(0.5*(self.size + (self.size/self.cells
                                                              if np.array(center).dtype in np.sctypes['int'] else
                                                              0)) if periodic else c))
        coords_rot = R.broadcast_to(tuple(self.cells))@coords

        with np.errstate(all='ignore'):
            mask = np.sum(np.power(coords_rot/r,2.0**np.array(exponent)),axis=-1) > 1.0

        if periodic:                                                                                # translate back to center
            mask = np.roll(mask,((c/self.size-0.5)*self.cells).round().astype(int),(0,1,2))

        return Grid(material = np.where(np.logical_not(mask) if inverse else mask,
                                        self.material,
                                        np.nanmax(self.material)+1 if fill is None else fill),
                    size     = self.size,
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Grid','add_primitive')],
                   )


    def mirror(self,directions,reflect=False):
        """
        Mirror grid along given directions.

        Parameters
        ----------
        directions : iterable containing str
            Direction(s) along which the grid is mirrored.
            Valid entries are 'x', 'y', 'z'.
        reflect : bool, optional
            Reflect (include) outermost layers. Defaults to False.

        """
        valid = ['x','y','z']
        if not set(directions).issubset(valid):
            raise ValueError(f'invalid direction {set(directions).difference(valid)} specified')

        limits = [None,None] if reflect else [-2,0]
        mat = self.material.copy()

        if 'x' in directions:
            mat = np.concatenate([mat,mat[limits[0]:limits[1]:-1,:,:]],0)
        if 'y' in directions:
            mat = np.concatenate([mat,mat[:,limits[0]:limits[1]:-1,:]],1)
        if 'z' in directions:
            mat = np.concatenate([mat,mat[:,:,limits[0]:limits[1]:-1]],2)

        return Grid(material = mat,
                    size     = self.size/self.cells*np.asarray(mat.shape),
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Grid','mirror')],
                   )


    def flip(self,directions):
        """
        Flip grid along given directions.

        Parameters
        ----------
        directions : iterable containing str
            Direction(s) along which the grid is flipped.
            Valid entries are 'x', 'y', 'z'.

        """
        valid = ['x','y','z']
        if not set(directions).issubset(valid):
            raise ValueError(f'invalid direction {set(directions).difference(valid)} specified')

        mat = np.flip(self.material, (valid.index(d) for d in directions if d in valid))

        return Grid(material = mat,
                    size     = self.size,
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Grid','flip')],
                   )


    def scale(self,cells,periodic=True):
        """
        Scale grid to new cells.

        Parameters
        ----------
        cells : numpy.ndarray of shape (3)
            Number of cells in x,y,z direction.
        periodic : Boolean, optional
            Assume grid to be periodic. Defaults to True.

        """
        return Grid(material = ndimage.interpolation.zoom(
                                                           self.material,
                                                           cells/self.cells,
                                                           output=self.material.dtype,
                                                           order=0,
                                                           mode=('wrap' if periodic else 'nearest'),
                                                           prefilter=False
                                                          ),
                    size     = self.size,
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Grid','scale')],
                   )


    def clean(self,stencil=3,selection=None,periodic=True):
        """
        Smooth grid by selecting most frequent material index within given stencil at each location.

        Parameters
        ----------
        stencil : int, optional
            Size of smoothing stencil.
        selection : list, optional
            Field values that can be altered. Defaults to all.
        periodic : Boolean, optional
            Assume grid to be periodic. Defaults to True.

        """
        def mostFrequent(arr,selection=None):
            me = arr[arr.size//2]
            if selection is None or me in selection:
                unique, inverse = np.unique(arr, return_inverse=True)
                return unique[np.argmax(np.bincount(inverse))]
            else:
                return me

        return Grid(material = ndimage.filters.generic_filter(
                                                               self.material,
                                                               mostFrequent,
                                                               size=(stencil if selection is None else stencil//2*2+1,)*3,
                                                               mode=('wrap' if periodic else 'nearest'),
                                                               extra_keywords=dict(selection=selection),
                                                              ).astype(self.material.dtype),
                    size     = self.size,
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Grid','clean')],
                   )


    def renumber(self):
        """Renumber sorted material indices as 0,...,N-1."""
        _,renumbered = np.unique(self.material,return_inverse=True)

        return Grid(material = renumbered.reshape(self.cells),
                    size     = self.size,
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Grid','renumber')],
                   )


    def rotate(self,R,fill=None):
        """
        Rotate grid (pad if required).

        Parameters
        ----------
        R : damask.Rotation
            Rotation to apply to the grid.
        fill : int or float, optional
            Material index to fill the corners. Defaults to material.max() + 1.

        """
        if fill is None: fill = np.nanmax(self.material) + 1
        dtype = float if isinstance(fill,float) or self.material.dtype in np.sctypes['float'] else int

        material = self.material
        # These rotations are always applied in the reference coordinate system, i.e. (z,x,z) not (z,x',z'')
        # see https://www.cs.utexas.edu/~theshark/courses/cs354/lectures/cs354-14.pdf
        for angle,axes in zip(R.as_Euler_angles(degrees=True)[::-1], [(0,1),(1,2),(0,1)]):
            material_temp = ndimage.rotate(material,angle,axes,order=0,prefilter=False,output=dtype,cval=fill)
            # avoid scipy interpolation errors for rotations close to multiples of 90°
            material = material_temp if np.prod(material_temp.shape) != np.prod(material.shape) else \
                       np.rot90(material,k=np.rint(angle/90.).astype(int),axes=axes)

        origin = self.origin-(np.asarray(material.shape)-self.cells)*.5 * self.size/self.cells

        return Grid(material = material,
                    size     = self.size/self.cells*np.asarray(material.shape),
                    origin   = origin,
                    comments = self.comments+[util.execution_stamp('Grid','rotate')],
                   )


    def canvas(self,cells=None,offset=None,fill=None):
        """
        Crop or enlarge/pad grid.

        Parameters
        ----------
        cells : numpy.ndarray of shape (3)
            Number of cells  x,y,z direction.
        offset : numpy.ndarray of shape (3)
            Offset (measured in cells) from old to new grid [0,0,0].
        fill : int or float, optional
            Material index to fill the background. Defaults to material.max() + 1.

        """
        if offset is None: offset = 0
        if fill is None: fill = np.nanmax(self.material) + 1
        dtype = float if int(fill) != fill or self.material.dtype in np.sctypes['float'] else int

        canvas = np.full(self.cells if cells is None else cells,fill,dtype)

        LL = np.clip( offset,           0,np.minimum(self.cells,     cells+offset))
        UR = np.clip( offset+cells,     0,np.minimum(self.cells,     cells+offset))
        ll = np.clip(-offset,           0,np.minimum(     cells,self.cells-offset))
        ur = np.clip(-offset+self.cells,0,np.minimum(     cells,self.cells-offset))

        canvas[ll[0]:ur[0],ll[1]:ur[1],ll[2]:ur[2]] = self.material[LL[0]:UR[0],LL[1]:UR[1],LL[2]:UR[2]]

        return Grid(material = canvas,
                    size     = self.size/self.cells*np.asarray(canvas.shape),
                    origin   = self.origin+offset*self.size/self.cells,
                    comments = self.comments+[util.execution_stamp('Grid','canvas')],
                   )


    def substitute(self,from_material,to_material):
        """
        Substitute material indices.

        Parameters
        ----------
        from_material : iterable of ints
            Material indices to be substituted.
        to_material : iterable of ints
            New material indices.

        """
        def mp(entry,mapper):
            return mapper[entry] if entry in mapper else entry

        mp = np.vectorize(mp)
        mapper = dict(zip(from_material,to_material))

        return Grid(material = mp(self.material,mapper).reshape(self.cells),
                    size     = self.size,
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Grid','substitute')],
                   )


    def sort(self):
        """Sort material indices such that min(material) is located at (0,0,0)."""
        a = self.material.flatten(order='F')
        from_ma = pd.unique(a)
        sort_idx = np.argsort(from_ma)
        ma = np.unique(a)[sort_idx][np.searchsorted(from_ma,a,sorter = sort_idx)]

        return Grid(material = ma.reshape(self.cells,order='F'),
                    size     = self.size,
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Grid','sort')],
                   )


    def vicinity_offset(self,vicinity=1,offset=None,trigger=[],periodic=True):
        """
        Offset material index of points in the vicinity of xxx.

        Different from themselves (or listed as triggers) within a given (cubic) vicinity,
        i.e. within the region close to a grain/phase boundary.
        ToDo: use include/exclude as in seeds.from_grid

        Parameters
        ----------
        vicinity : int, optional
            Voxel distance checked for presence of other materials.
            Defaults to 1.
        offset : int, optional
            Offset (positive or negative) to tag material indices,
            defaults to material.max()+1.
        trigger : list of ints, optional
            List of material indices that trigger a change.
            Defaults to [], meaning that any different neighbor triggers a change.
        periodic : Boolean, optional
            Assume grid to be periodic. Defaults to True.

        """
        def tainted_neighborhood(stencil,trigger):

            me = stencil[stencil.shape[0]//2]
            return np.any(stencil != me if len(trigger) == 0 else
                          np.in1d(stencil,np.array(list(set(trigger) - {me}))))

        offset_ = np.nanmax(self.material)+1 if offset is None else offset
        mask = ndimage.filters.generic_filter(self.material,
                                              tainted_neighborhood,
                                              size=1+2*vicinity,
                                              mode='wrap' if periodic else 'nearest',
                                              extra_keywords={'trigger':trigger})

        return Grid(material = np.where(mask, self.material + offset_,self.material),
                    size     = self.size,
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Grid','vicinity_offset')],
                   )


    def get_grain_boundaries(self,periodic=True,directions='xyz'):
        """
        Create VTK unstructured grid containing grain boundaries.

        Parameters
        ----------
        periodic : Boolean, optional
            Assume grid to be periodic. Defaults to True.
        directions : iterable containing str, optional
            Direction(s) along which the boundaries are determined.
            Valid entries are 'x', 'y', 'z'. Defaults to 'xyz'.

        """
        valid = ['x','y','z']
        if not set(directions).issubset(valid):
            raise ValueError(f'invalid direction {set(directions).difference(valid)} specified')

        o = [[0, self.cells[0]+1,           np.prod(self.cells[:2]+1)+self.cells[0]+1, np.prod(self.cells[:2]+1)],
             [0, np.prod(self.cells[:2]+1), np.prod(self.cells[:2]+1)+1,               1],
             [0, 1,                         self.cells[0]+1+1,                         self.cells[0]+1]] # offset for connectivity

        connectivity = []
        for i,d in enumerate(['x','y','z']):
            if d not in directions: continue
            mask = self.material != np.roll(self.material,1,i)
            for j in [0,1,2]:
                mask = np.concatenate((mask,np.take(mask,[0],j)*(i==j)),j)
            if i == 0 and not periodic: mask[0,:,:] = mask[-1,:,:] = False
            if i == 1 and not periodic: mask[:,0,:] = mask[:,-1,:] = False
            if i == 2 and not periodic: mask[:,:,0] = mask[:,:,-1] = False

            base_nodes = np.argwhere(mask.flatten(order='F')).reshape(-1,1)
            connectivity.append(np.block([base_nodes + o[i][k] for k in range(4)]))

        coords = grid_filters.coordinates0_node(self.cells,self.size,self.origin).reshape(-1,3,order='F')
        return VTK.from_unstructured_grid(coords,np.vstack(connectivity),'QUAD')
