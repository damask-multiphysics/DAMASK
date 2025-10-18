import os
import copy
import numbers
import multiprocessing as mp
from functools import partial
from typing import Optional, Union, Sequence
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
from scipy import ndimage, spatial, interpolate

from . import VTK
from . import util
from . import grid_filters
from . import Rotation
from . import Table
from . import Colormap
from ._typehints import FloatSequence, IntSequence, NumpyRngSeed


class IcDict(dict):
    """Dict that validates and broadcasts initial conditions on assignment."""

    def __init__(self, cells: tuple[int, int, int]):
        """
        New initial condition dictionary.

        Parameters
        ----------
        cells : tuple of int, len (3)
            Cell counts along x,y,z direction.
        """
        super().__init__()
        self.cells = cells

    def __setitem__(self, k: str, v: Union[int, float, np.ndarray]):
        """
        Set a single initial condition key with validation and broadcasting.

        Parameters
        ----------
        k : str
            Name of the initial condition field.
        v : int, float, or np.ndarray
            Value to assign to the field. Allowed types and shapes:
            - Scalars (`int` or `float`)
              → broadcast to `cells`.
            - Small arrays (`np.ndarray`) of shape (), (1,), or (3,)
              → broadcast to `cells + value.shape`.
            - Full-field arrays (`np.ndarray`) with leading dimensions matching `cells` and trailing shape (), (1,), or (3,)
              → copied as-is.

        Raises
        ------
        ValueError
            If `v` has an invalid shape that does not match the allowed rules.

        Notes
        -----
        Broadcasting ensures that all initial condition fields conform to the
        grid defined by `cells`.
        This method is called both for full-property assignment (via the parent setter)
        and per-key assignment (`GeomGrid.initial_conditions[key] = value`), so all
        validation and broadcasting logic is centralized here.
        """
        if not (    isinstance(v, numbers.Real)
                or (isinstance(v, np.ndarray) and v.shape in [(), (1,), (3,)])
                or (isinstance(v, np.ndarray) and len(v.shape) >= 3 and
                    v.shape[:3] == self.cells and v.shape[3:] in [(), (1,), (3,)])
               ):
            raise ValueError(f'initial condition "{k}" must be [a field of] scalars or three-dimensional vectors')

        super().__setitem__(k,
                            v if isinstance(v, np.ndarray) and len(v.shape) >= 3 and v.shape[:3] == self.cells else
                            np.broadcast_to(v, self.cells + v.shape) if isinstance(v, np.ndarray) else
                            np.broadcast_to(v, self.cells)
                            )


class GeomGrid:
    """
    Geometry definition for grid solvers.

    Create and manipulate geometry definitions for storage as VTK ImageData
    files ('.vti' extension). A grid has a physical size, a coordinate origin,
    and contains the material ID (indexing an entry in 'material.yaml')
    as well as initial condition fields.
    """

    def __init__(self,
                 material: np.ndarray,
                 size: FloatSequence,
                 origin: FloatSequence = np.zeros(3),
                 initial_conditions: Optional[dict[str,np.ndarray]] = None,
                 comments: Union[None, str, Sequence[str]] = None):
        """
        New geometry definition for grid solvers.

        Parameters
        ----------
        material : numpy.ndarray of int, shape (:,:,:)
            Material indices. The shape of the material array defines
            the number of cells.
        size : sequence of float, len (3)
            Physical size of grid in meter.
        origin : sequence of float, len (3), optional
            Coordinates of grid origin in meter. Defaults to [0.0,0.0,0.0].
        initial_conditions : dictionary, optional
            Initial condition label and field values at each grid point.
            If leading shape deviates from material shape,
            the (constant) value is broadcast across the grid.
        comments : (sequence of) str, optional
            Additional, human-readable information, e.g. history of operations.
        """
        self.material = material
        self.size = size
        self.origin = origin
        self._ic = IcDict(tuple(self.cells))
        self.initial_conditions = {} if initial_conditions is None else initial_conditions
        self.comments = [] if comments is None else \
                        [comments] if isinstance(comments,str) else \
                        [str(c) for c in comments]

    def __repr__(self) -> str:
        """
        Return repr(self).

        Give short, human-readable summary.
        """
        mat_min = np.nanmin(self.material)
        mat_max = np.nanmax(self.material)
        mat_N   = self.N_materials
        return util.srepr([
               f'cells:  {util.srepr(self.cells, " × ")}',
               f'size:   {util.srepr(self.size,  " × ")} m³',
               f'origin: {util.srepr(self.origin,"   ")} m',
               f'# materials: {mat_N}' + ('' if mat_min == 0 and mat_max == mat_N-1 else
                                          f' (min: {mat_min}, max: {mat_max})')
               ]+(['initial_conditions:']+[f'  - {f}'+(f' {data.shape[3:]}' if len(data.shape)>3 else '')
                                           for f,data in self.initial_conditions.items()] if self.initial_conditions else []))


    def __copy__(self) -> 'GeomGrid':
        """
        Return deepcopy(self).

        Create deep copy.
        """
        return copy.deepcopy(self)

    copy = __copy__


    def __eq__(self,
               other: object) -> bool:
        """
        Return self==other.

        Test equality of other.

        Parameters
        ----------
        other : damask.GeomGrid
            GeomGrid to compare self against.

        Returns
        -------
        equal : bool
            Whether both arguments are equal.
        """
        if not isinstance(other, GeomGrid):
            return NotImplemented
        return self.match(other) and bool(    other.initial_conditions.keys() == self.initial_conditions.keys()
                                          and all(np.allclose(other.initial_conditions[k], self.initial_conditions[k])
                                                  for k in self.initial_conditions)
                                         )


    def match(self,
              other: object) -> bool:
        """
        Test geometric equality of other, i.e., ignoring initial conditions.

        Parameters
        ----------
        other : damask.GeomGrid
            GeomGrid to compare self against.

        Returns
        -------
        match : bool
            Whether both arguments are geometrically equal.

        Notes
        -----
        This comparison does not consider initial conditions.
        """
        if not isinstance(other, GeomGrid):
            return NotImplemented
        return bool(    np.allclose(other.size,self.size)
                    and np.allclose(other.origin,self.origin)
                    and np.all(other.cells == self.cells)
                    and np.all(other.material == self.material)
                   )


    @property
    def material(self) -> np.ndarray:
        """Material indices."""
        return self._material

    @material.setter
    def material(self,
                 material: np.ndarray):
        if len(material.shape) != 3:
            raise ValueError(f'invalid material shape {material.shape}')
        if material.dtype not in [np.float32,np.float64, np.int32,np.int64]:
            raise TypeError(f'invalid material data type "{material.dtype}"')

        self._material = np.copy(material)

        if self.material.dtype in [np.float32,np.float64] and \
           np.all(self.material == self.material.astype(np.int64).astype(float)):
            self._material = self.material.astype(np.int64)


    @property
    def size(self) -> np.ndarray:
        """Edge lengths of grid in meter."""
        return self._size

    @size.setter
    def size(self,
             size: FloatSequence):
        if len(size) != 3 or any(np.asarray(size) < 0):
            raise ValueError(f'invalid size {size}')

        self._size = np.array(size,np.float64)

    @property
    def origin(self) -> np.ndarray:
        """Vector to grid origin in meter."""
        return self._origin

    @origin.setter
    def origin(self,
               origin: FloatSequence):
        if len(origin) != 3:
            raise ValueError(f'invalid origin {origin}')

        self._origin = np.array(origin,np.float64)


    @property
    def initial_conditions(self) -> dict[str, np.ndarray]:
        """Fields of initial conditions."""
        return self._ic

    @initial_conditions.setter
    def initial_conditions(self, ic: dict[str, np.ndarray]):
        """Set initial conditions, broadcasting to the cell grid shape if necessary."""
        if not isinstance(ic, dict):
            raise TypeError('initial conditions must be a dictionary')

        self._ic.clear()
        for k, v in ic.items():
            self._ic[k] = v


    @property
    def cells(self) -> np.ndarray:
        """Cell counts along x,y,z direction."""
        return np.asarray(self.material.shape)


    @property
    def N_materials(self) -> int:
        """Number of (unique) material indices within grid."""
        return np.unique(self.material).size


    @staticmethod
    def _load(fname: Union[str, Path], label: str) -> 'GeomGrid':
        """
        Load from VTK ImageData file.

        Parameters
        ----------
        fname : str or pathlib.Path
            VTK ImageData file to read.
            Valid extension is .vti, which will be appended if not given.
        label : str
            Label of the dataset containing the material IDs.

        Returns
        -------
        loaded : damask.GeomGrid
            Grid-based geometry from file.
        """
        v = VTK.load(fname if str(fname).endswith('.vti') else str(fname)+'.vti')
        cells = tuple(np.array(v.vtk_data.GetDimensions())-1)                                       # type: ignore[attr-defined]
        bbox  = np.array(v.vtk_data.GetBounds()).reshape(3,2).T
        ic = {l:v.get(l).reshape(cells+v.get(l).shape[1:],order='F') for l in set(v.labels['Cell Data']) - {label}}

        return GeomGrid(material = v.get(label).reshape(cells,order='F'),
                        size     = bbox[1] - bbox[0],
                        origin   = bbox[0],
                        initial_conditions = ic,
                        comments = v.comments,
                       )

    @staticmethod
    def load(fname: Union[str, Path]) -> 'GeomGrid':
        """
        Load from VTK ImageData file with material IDs stored as 'material'.

        Parameters
        ----------
        fname : str or pathlib.Path
            GeomGrid file to read.
            Valid extension is .vti, which will be appended if not given.

        Returns
        -------
        loaded : damask.GeomGrid
            Grid-based geometry from file.
        """
        return GeomGrid._load(fname,'material')


    @staticmethod
    def load_SPPARKS(fname: Union[str, Path]) -> 'GeomGrid':
        """
        Load from SPPARKS VTK dump.

        Parameters
        ----------
        fname : str or pathlib.Path
            SPPARKS VTK dump file to read.
            Valid extension is .vti, which will be appended if not given.

        Returns
        -------
        loaded : damask.GeomGrid
            Grid-based geometry from file.

        Notes
        -----
        A SPPARKS VTI dump is equivalent to a DAMASK VTI file,
        but stores the materialID information as 'Spin' rather than 'material'.
        """
        return GeomGrid._load(fname,'Spin')


    @staticmethod
    def load_Neper(fname: Union[str, Path]) -> 'GeomGrid':
        """
        Load from Neper VTK file.

        Parameters
        ----------
        fname : str or pathlib.Path
            Geometry file to read.

        Returns
        -------
        loaded : damask.GeomGrid
            Grid-based geometry from file.

        Notes
        -----
        Material indices in Neper usually start at 1 unless
        a buffer material with index 0 is added.

        Examples
        --------
        Read a periodic polycrystal generated with Neper.

        >>> import damask
        >>> N_grains = 20
        >>> cells = (32,32,32)
        >>> damask.util.run(f'neper -T -n {N_grains} -tesrsize {cells[0]}:{cells[1]}:{cells[2]} -periodicity all -format vtk')
        stdioTuple(stdout=...
        >>> damask.GeomGrid.load_Neper(f'n{N_grains}-id1.vtk').renumber()
        cells:  32 × 32 × 32
        size:   1.0 × 1.0 × 1.0 m³
        origin: 0.0   0.0   0.0 m
        # materials: 20
        """
        v = VTK.load(fname,'ImageData')
        cells = np.array(v.vtk_data.GetDimensions())-1                                              # type: ignore[attr-defined]
        bbox  = np.array(v.vtk_data.GetBounds()).reshape(3,2).T

        return GeomGrid(material = v.get('MaterialId').reshape(cells,order='F').astype('int32',casting='unsafe'),
                        size     = bbox[1] - bbox[0],
                        origin   = bbox[0],
                        comments = util.execution_stamp('GeomGrid','load_Neper'),
                       )


    @staticmethod
    def load_DREAM3D(fname: Union[str, Path],
                     feature_IDs: Optional[str] = None,
                     cell_data: Optional[str] = None,
                     phases: str = 'Phases',
                     Euler_angles: str = 'EulerAngles',
                     base_group: Optional[str] = None) -> 'GeomGrid':
        """
        Load DREAM.3D (HDF5) file.

        Data in DREAM.3D files can be stored per cell ('CellData')
        and/or per grain ('Grain Data'). Per default, i.e. if
        'feature_IDs' is None, cell-wise data is assumed.

        Parameters
        ----------
        fname : str or pathlib.Path
            Filename of the DREAM.3D (HDF5) file.
        feature_IDs : str, optional
            Name of the dataset containing the mapping between cells and
            grain-wise data. Defaults to 'None', in which case cell-wise
            data is used.
        cell_data : str, optional
            Name of the group (folder) containing cell-wise data. Defaults to
            None in wich case it is automatically detected.
        phases : str, optional
            Name of the dataset containing the phase ID. It is not used for
            grain-wise data, i.e. when feature_IDs is not None.
            Defaults to 'Phases'.
        Euler_angles : str, optional
            Name of the dataset containing the crystallographic orientation as
            Euler angles in radians It is not used for grain-wise data, i.e.
            when feature_IDs is not None. Defaults to 'EulerAngles'.
        base_group : str, optional
            Path to the group (folder) that contains geometry (_SIMPL_GEOMETRY),
            and grain- or cell-wise data. Defaults to None, in which case
            it is set as the path that contains _SIMPL_GEOMETRY/SPACING.

        Returns
        -------
        loaded : damask.GeomGrid
            Grid-based geometry from file.

        Notes
        -----
        A grain-wise geometry definition is based on segmented data from the
        DREAM.3D file. This data is typically available when the microstructure
        was synthetically created. In cell-wise representations, cells having
        the same orientation and phase are grouped. Since synthetically created
        microstructures have typically no in-grain scatter, cell-wise grids
        can appear to be segmented.

        damask.ConfigMaterial.load_DREAM3D creates the corresponding
        material definition. Since the numbering of materials in cell-wise
        and grain-wise grids is different, it is imperative to use the same
        mode for both load_DREAM3D functions. That means, if the "feature_IDs"
        argument is used for this function, the correct material configuration
        is only obtained if the "grain_data" argument is used when calling
        damask.ConfigMaterial.load_DREAM3D.

        Versions 8.0 and later of the DREAM.3D file format are not yet supported.
        """
        with h5py.File(fname, 'r') as f:
            if (file_version := util.version(f.attrs['FileVersion'].decode()+'.0')) > '7.0.0':
                raise ValueError(f'DREAM.3D file format {file_version} is not supported')
            b = util.DREAM3D_base_group(f) if base_group is None else base_group
            c = util.DREAM3D_cell_data_group(f) if cell_data is None else cell_data

            cells  = f['/'.join([b,'_SIMPL_GEOMETRY','DIMENSIONS'])][()]
            size   = f['/'.join([b,'_SIMPL_GEOMETRY','SPACING'])] * cells
            origin = f['/'.join([b,'_SIMPL_GEOMETRY','ORIGIN'])][()]

            if feature_IDs is None:
                phase = f['/'.join([b,c,phases])][()].reshape(-1,1)
                O = Rotation.from_Euler_angles(f['/'.join([b,c,Euler_angles])]).as_quaternion().reshape(-1,4) # noqa
                unique,unique_inverse = np.unique(np.hstack([O,phase]),return_inverse=True,axis=0)
                ma = np.arange(cells.prod()) if len(unique) == cells.prod() else \
                     np.arange(unique.size)[np.argsort(pd.unique(unique_inverse.squeeze()))][unique_inverse]
            else:
                ma = f['/'.join([b,c,feature_IDs])][()].flatten()

        return GeomGrid(material = ma.reshape(cells,order='F'),
                        size     = size,
                        origin   = origin,
                        comments = util.execution_stamp('GeomGrid','load_DREAM3D'),
                       )


    @staticmethod
    def from_table(table: Table,
                   coordinates: str,
                   labels: Union[str, Sequence[str]],
                   atol: float = 0.0) -> 'GeomGrid':
        """
        Create grid from ASCII table.

        Parameters
        ----------
        table : damask.Table
            Table that contains material information.
        coordinates : str
            Label of the vector column containing the spatial coordinates.
            Need to be ordered (1./x fast, 3./z slow).
        labels : (sequence of) str
            Label(s) of the columns containing the material definition.
            Each unique combination of values results in one material ID.
        atol : float, optional
            Absolute tolerance to consider grid coordinates equivalent.
            Defaults to 0.0.

        Returns
        -------
        new : damask.GeomGrid
            Grid-based geometry from values in table.
        """
        cells,size,origin = grid_filters.cellsSizeOrigin_coordinates0_point(table.get(coordinates),atol=atol)

        unique,inverse = table[labels].unique(return_inverse=True)

        return GeomGrid(material = np.arange(len(unique))[inverse].reshape(cells,order='F'),
                        size     = size,
                        origin   = origin,
                        comments = util.execution_stamp('GeomGrid','from_table'),
                       )


    @staticmethod
    def _find_closest_seed(seeds: np.ndarray,
                           weights: np.ndarray,
                           point: np.ndarray) -> np.integer:
        return np.argmin(np.sum((np.broadcast_to(point,(len(seeds),3))-seeds)**2,axis=1) - weights)

    @staticmethod
    def from_Laguerre_tessellation(cells: IntSequence,
                                   size: FloatSequence,
                                   seeds: np.ndarray,
                                   weights: FloatSequence,
                                   material: Optional[IntSequence] = None,
                                   periodic: bool = True):
        """
        Create grid from Laguerre tessellation.

        Parameters
        ----------
        cells : sequence of int, len (3)
            Cell counts along x,y,z direction.
        size : sequence of float, len (3)
            Edge lengths of the grid in meter.
        seeds : numpy.ndarray of float, shape (:,3)
            Position of the seed points in meter. All points need
            to lay within the box [(0,0,0),size].
        weights : sequence of float, len (seeds.shape[0])
            Weights of the seeds. Setting all weights to 1.0 gives a
            standard Voronoi tessellation.
        material : sequence of int, len (seeds.shape[0]), optional
            Material ID of the seeds.
            Defaults to None, in which case materials are consecutively numbered.
        periodic : bool, optional
            Assume grid to be periodic. Defaults to True.

        Returns
        -------
        new : damask.GeomGrid
            Grid-based geometry from tessellation.

        Notes
        -----
        damask.seeds contains functionality for seed generation.
        """
        weights_p: FloatSequence
        if periodic:
            weights_p = np.tile(weights,27)                                                         # Laguerre weights (1,2,3,1,2,3,...,1,2,3)
            seeds_p = np.vstack((seeds  -np.array([size[0],0.,0.]),seeds,  seeds  +np.array([size[0],0.,0.])))
            seeds_p = np.vstack((seeds_p-np.array([0.,size[1],0.]),seeds_p,seeds_p+np.array([0.,size[1],0.])))
            seeds_p = np.vstack((seeds_p-np.array([0.,0.,size[2]]),seeds_p,seeds_p+np.array([0.,0.,size[2]])))
        else:
            weights_p = np.array(weights,float)
            seeds_p   = seeds

        coords = grid_filters.coordinates0_point(cells,size).reshape(-1,3)

        pool = mp.Pool(int(os.environ.get('OMP_NUM_THREADS',4)))
        result = pool.map_async(partial(GeomGrid._find_closest_seed,seeds_p,weights_p), coords)
        pool.close()
        pool.join()
        material_ = np.array(result.get()).reshape(cells)

        if periodic: material_ %= len(weights)

        return GeomGrid(material = material_ if material is None else np.array(material)[material_],
                        size     = size,
                        comments = util.execution_stamp('GeomGrid','from_Laguerre_tessellation'),
                       )


    @staticmethod
    def from_Voronoi_tessellation(cells: IntSequence,
                                  size: FloatSequence,
                                  seeds: np.ndarray,
                                  material: Optional[IntSequence] = None,
                                  periodic: bool = True) -> 'GeomGrid':
        """
        Create grid from Voronoi tessellation.

        Parameters
        ----------
        cells : sequence of int, len (3)
            Cell counts along x,y,z direction.
        size : sequence of float, len (3)
            Edge lengths of the grid in meter.
        seeds : numpy.ndarray of float, shape (:,3)
            Position of the seed points in meter. All points need
            to lay within the box [(0,0,0),size].
        material : sequence of int, len (seeds.shape[0]), optional
            Material ID of the seeds.
            Defaults to None, in which case materials are consecutively numbered.
        periodic : bool, optional
            Assume grid to be periodic. Defaults to True.

        Returns
        -------
        new : damask.GeomGrid
            Grid-based geometry from tessellation.

        Notes
        -----
        damask.seeds contains functionality for seed generation.

        Examples
        --------
        Generate microstructure with three grains.

        >>> import numpy as np
        >>> import damask
        >>> seeds = np.array([[0.30828, 0.64020, 0.65237],
        ...                   [0.62200, 0.56858, 0.32842],
        ...                   [0.57315, 0.94534, 0.87531]])*1e-6
        >>> damask.GeomGrid.from_Voronoi_tessellation(cells=[10,10,10],
        ...                                           size=[1e-6]*3,
        ...                                           seeds=seeds)
        cells:  10 × 10 × 10
        size:   1e-06 × 1e-06 × 1e-06 m³
        origin: 0.0   0.0   0.0 m
        # materials: 3
        """
        coords = grid_filters.coordinates0_point(cells,size).reshape(-1,3)
        tree = spatial.KDTree(seeds,boxsize=np.asarray(size) if periodic else None)
        material_ = tree.query(coords, workers = int(os.environ.get('OMP_NUM_THREADS',4)))[1]

        return GeomGrid(material = (material_ if material is None else np.array(material)[material_]).reshape(cells),
                        size     = size,
                        comments = util.execution_stamp('GeomGrid','from_Voronoi_tessellation'),
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
    def from_minimal_surface(cells: IntSequence,
                             size: FloatSequence,
                             surface: str,
                             threshold: float = 0.0,
                             periods: int = 1,
                             materials: IntSequence = (0,1)) -> 'GeomGrid':
        """
        Create grid from definition of triply-periodic minimal surface.

        Parameters
        ----------
        cells : sequence of int, len (3)
            Cell counts along x,y,z direction.
        size : sequence of float, len (3)
            Edge lengths of the grid in meter.
        surface : str
            Type of the minimal surface. See notes for details.
        threshold : float, optional
            Threshold of the minimal surface. Defaults to 0.0.
        periods : int, optional
            Number of periods per unit cell. Defaults to 1.
        materials : sequence of int, len (2)
            Material IDs. Defaults to (0,1).

        Returns
        -------
        new : damask.GeomGrid
            Grid-based geometry from definition of minimal surface.

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

        Examples
        --------
        Minimal surface of 'Gyroid' type.

        >>> import numpy as np
        >>> import damask
        >>> damask.GeomGrid.from_minimal_surface([64]*3,np.ones(3)*1.e-4,'Gyroid')
        cells:  64 × 64 × 64
        size:   0.0001 × 0.0001 × 0.0001 m³
        origin: 0.0   0.0   0.0 m
        # materials: 2

        Minimal surface of 'Neovius' type with specific material IDs.

        >>> import numpy as np
        >>> import damask
        >>> damask.GeomGrid.from_minimal_surface([80]*3,np.ones(3)*5.e-4,
        ...                                  'Neovius',materials=(1,5))
        cells:  80 × 80 × 80
        size:   0.0005 × 0.0005 × 0.0005 m³
        origin: 0.0   0.0   0.0 m
        # materials: 2 (min: 1, max: 5)
        """
        x,y,z = np.meshgrid(periods*2.0*np.pi*(np.arange(cells[0])+0.5)/cells[0],
                            periods*2.0*np.pi*(np.arange(cells[1])+0.5)/cells[1],
                            periods*2.0*np.pi*(np.arange(cells[2])+0.5)/cells[2],
                            indexing='ij',sparse=True)
        return GeomGrid(material = np.where(threshold < GeomGrid._minimal_surface[surface](x,y,z),materials[1],materials[0]),
                        size     = size,
                        comments = util.execution_stamp('GeomGrid','from_minimal_surface'),
                       )


    def save(self,
             fname: Union[str, Path],
             compress: bool = True):
        """
        Save as VTK ImageData file.

        Parameters
        ----------
        fname : str or pathlib.Path
            Filename to write.
            Valid extension is .vti, which will be appended if not given.
        compress : bool, optional
            Compress with zlib algorithm. Defaults to True.
        """
        v = VTK.from_image_data(self.cells,self.size,self.origin)\
               .set('material',self.material.flatten(order='F'))
        for label,data in self.initial_conditions.items():
            v = v.set(label,data.reshape((-1,)+data.shape[3:],order='F'))
        v.comments = self.comments

        v.save(fname,parallel=False,compress=compress)



    def show(self,
             colormap: Union[Colormap, str] = 'cividis') -> None:
        """
        Show on screen.

        Parameters
        ----------
        colormap : damask.Colormap or str, optional
            Colormap for visualization of material IDs. Defaults to 'cividis'.
        """
        VTK.from_image_data(self.cells,self.size,self.origin) \
           .set('material',self.material.flatten('F'),) \
           .show('material',colormap)


    def canvas(self,
               cells: Optional[IntSequence] = None,
               offset: Optional[IntSequence] = None,
               fill: Optional[int] = None) -> 'GeomGrid':
        """
        Crop or enlarge/pad grid.

        Parameters
        ----------
        cells : sequence of int, len (3), optional
            Cell counts along x,y,z direction.
        offset : sequence of int, len (3), optional
            Offset (measured in cells) from old to new grid.
            Defaults to [0,0,0].
        fill : int, optional
            Material ID to fill the background.
            Defaults to material.max()+1.

        Returns
        -------
        updated : damask.GeomGrid
            Updated grid-based geometry.

        Notes
        -----
        Existing initial condition fields will be removed.

        Examples
        --------
        Remove lower 1/2 of the microstructure in z-direction.

        >>> import numpy as np
        >>> import damask
        >>> g = damask.GeomGrid(np.zeros([32]*3,int),np.ones(3)*1e-3)
        >>> g.canvas([32,32,16],[0,0,16])
        cells:  32 × 32 × 16
        size:   0.001 × 0.001 × 0.0005 m³
        origin: 0.0   0.0   0.0005 m
        # materials: 1
        """
        offset_ = np.array(offset,np.int64) if offset is not None else np.zeros(3,np.int64)
        cells_ = np.array(cells,np.int64) if cells is not None else self.cells

        canvas = np.full(cells_,np.nanmax(self.material) + 1 if fill is None else fill,self.material.dtype)

        LL = np.clip( offset_,           0,np.minimum(self.cells,     cells_+offset_))
        UR = np.clip( offset_+cells_,    0,np.minimum(self.cells,     cells_+offset_))
        ll = np.clip(-offset_,           0,np.minimum(     cells_,self.cells-offset_))
        ur = np.clip(-offset_+self.cells,0,np.minimum(     cells_,self.cells-offset_))

        canvas[ll[0]:ur[0],ll[1]:ur[1],ll[2]:ur[2]] = self.material[LL[0]:UR[0],LL[1]:UR[1],LL[2]:UR[2]]

        return GeomGrid(material = canvas,
                        size     = self.size/self.cells*np.asarray(canvas.shape),
                        origin   = self.origin+offset_*self.size/self.cells,
                        comments = self.comments+[util.execution_stamp('GeomGrid','canvas')],
                       )


    def mirror(self,
               directions: Sequence[str],
               reflect: bool = False) -> 'GeomGrid':
        """
        Mirror grid along given directions.

        Parameters
        ----------
        directions : (sequence of) {'x', 'y', 'z'}
            Direction(s) along which the grid is mirrored.
        reflect : bool, optional
            Reflect (include) outermost layers. Defaults to False.

        Returns
        -------
        updated : damask.GeomGrid
            Updated grid-based geometry.

        Examples
        --------
        Mirror along y-direction.

        >>> import numpy as np
        >>> import damask
        >>> (g := damask.GeomGrid(np.arange(4*5*6).reshape([4,5,6]),np.ones(3)))
        cells:  4 × 5 × 6
        size:   1.0 × 1.0 × 1.0 m³
        origin: 0.0   0.0   0.0 m
        # materials: 120
        >>> g.mirror('y')
        cells:  4 × 8 × 6
        size:   1.0 × 1.6 × 1.0 m³
        origin: 0.0   0.0   0.0 m
        # materials: 120

        Reflect along x- and y-direction.

        >>> g.mirror('xy',reflect=True)
        cells:  8 × 10 × 6
        size:   2.0 × 2.0 × 1.0 m³
        origin: 0.0   0.0   0.0 m
        # materials: 120

        Independence of mirroring order.

        >>> g.mirror('xy') == g.mirror(['y','x'])
        True
        """
        valid = 'xyz'
        if not set(directions).issubset(valid):
            raise ValueError(f'invalid direction "{set(directions).difference(valid)}" specified')

        limits: Sequence[Optional[int]] = [None,None] if reflect else [-2,0]
        selection = (slice(limits[0],limits[1],-1),slice(None),slice(None))
        mat = self.material.copy()
        ic = self.initial_conditions.copy()

        for i,d in enumerate(valid):
            if d in directions:
                mat = np.concatenate([mat,
                                      mat[selection[0-i],selection[1-i],selection[2-i]]],
                                      axis=i)
                for label in ic:
                    ic[label] = np.concatenate([ic[label],
                                                ic[label][selection[0-i],selection[1-i],selection[2-i]]],
                                                axis=i)

        return GeomGrid(material = mat,
                        size     = self.size/self.cells*np.asarray(mat.shape),
                        origin   = self.origin,
                        initial_conditions = ic,
                        comments = self.comments+[util.execution_stamp('GeomGrid','mirror')],
                       )


    def flip(self,
             directions: Sequence[str]) -> 'GeomGrid':
        """
        Flip grid along given directions.

        Parameters
        ----------
        directions : (sequence of) {'x', 'y', 'z'}
            Direction(s) along which the grid is flipped.

        Returns
        -------
        updated : damask.GeomGrid
            Updated grid-based geometry.

        Examples
        --------
        Invariance of flipping order.

        >>> import numpy as np
        >>> import damask
        >>> (g := damask.GeomGrid(np.arange(4*5*6).reshape([4,5,6]),np.ones(3)))
        cells:  4 × 5 × 6
        size:   1.0 × 1.0 × 1.0 m³
        origin: 0.0   0.0   0.0 m
        # materials: 120
        >>> g.flip('xyz') == g.flip(['x','z','y'])
        True

        Invariance of flipping a (fully) mirrored grid.

        >>> g.mirror('x',reflect=True) == g.mirror('x',reflect=True).flip('x')
        True
        """
        valid = 'xyz'
        if not set(directions).issubset(valid):
            raise ValueError(f'invalid direction "{set(directions).difference(valid)}" specified')

        mat = np.flip(self.material, [valid.index(d) for d in directions if d in valid])
        ic = {}
        for label in self.initial_conditions:
            ic[label] = np.flip(self.initial_conditions[label],[valid.index(d) for d in directions if d in valid])
        return GeomGrid(material = mat,
                        size     = self.size,
                        origin   = self.origin,
                        initial_conditions = ic,
                        comments = self.comments+[util.execution_stamp('GeomGrid','flip')],
                       )


    def rotate(self,
               R: Rotation,
               fill: Optional[int] = None) -> 'GeomGrid':
        """
        Rotate grid (possibly extending its bounding box).

        Parameters
        ----------
        R : damask.Rotation
            Rotation to apply to the grid.
        fill : int, optional
            Material ID to fill enlarged bounding box.
            Defaults to material.max()+1.

        Returns
        -------
        updated : damask.GeomGrid
            Updated grid-based geometry.

        Notes
        -----
        Existing initial condition fields will be removed.

        Examples
        --------
        Rotation by 180° (π) is equivalent to twice flipping.

        >>> import numpy as np
        >>> import damask
        >>> (g := damask.GeomGrid(np.arange(4*5*6).reshape([4,5,6]),np.ones(3)))
        cells:  4 × 5 × 6
        size:   1.0 × 1.0 × 1.0 m³
        origin: 0.0   0.0   0.0 m
        # materials: 120
        >>> g.rotate(damask.Rotation.from_axis_angle([0,0,1,180],degrees=True)) == g.flip('xy')
        True
        """
        material = self.material
        # These rotations are always applied in the reference coordinate system, i.e. (z,x,z) not (z,x',z'')
        # see https://www.cs.utexas.edu/~theshark/courses/cs354/lectures/cs354-14.pdf
        for angle,axes in zip(R.as_Euler_angles(degrees=True)[::-1], [(0,1),(1,2),(0,1)]):
            material_temp = ndimage.rotate(material,angle,axes,order=0,prefilter=False,
                                           output=self.material.dtype,
                                           cval=np.nanmax(self.material) + 1 if fill is None else fill)
            # avoid scipy interpolation errors for rotations close to multiples of 90°
            material = material_temp if np.prod(material_temp.shape) != np.prod(material.shape) else \
                       np.rot90(material,k=np.rint(angle/90.).astype(np.int64),axes=axes)

        origin = self.origin-(np.asarray(material.shape)-self.cells)*.5 * self.size/self.cells

        return GeomGrid(material = material,
                        size     = self.size/self.cells*np.asarray(material.shape),
                        origin   = origin,
                        comments = self.comments+[util.execution_stamp('GeomGrid','rotate')],
                       )


    def scale(self,
              cells: IntSequence) -> 'GeomGrid':
        """
        Scale grid to new cell counts.

        Parameters
        ----------
        cells : sequence of int, len (3)
            Cell counts along x,y,z direction.

        Returns
        -------
        updated : damask.GeomGrid
            Updated grid-based geometry.

        Examples
        --------
        Double grid resolution.

        >>> import numpy as np
        >>> import damask
        >>> (g := damask.GeomGrid(np.zeros([32]*3,int),np.ones(3)*1e-4))
        cells:  32 × 32 × 32
        size:   0.0001 × 0.0001 × 0.0001 m³
        origin: 0.0   0.0   0.0 m
        # materials: 1
        >>> g.scale(g.cells*2)
        cells:  64 × 64 × 64
        size:   0.0001 × 0.0001 × 0.0001 m³
        origin: 0.0   0.0   0.0 m
        # materials: 1
        """
        orig = tuple(map(np.linspace,self.origin             + self.size/self.cells*.5,
                                     self.origin + self.size - self.size/self.cells*.5,self.cells))
        interpolator = partial(interpolate.RegularGridInterpolator,
                               points=orig,method='nearest',bounds_error=False,fill_value=None)
        new = grid_filters.coordinates0_point(cells,self.size,self.origin)

        return GeomGrid(material = interpolator(values=self.material)(new).astype(int),
                        size     = self.size,
                        origin   = self.origin,
                        initial_conditions = {k: interpolator(values=v)(new)
                                              for k,v in self.initial_conditions.items()},
                        comments = self.comments+[util.execution_stamp('GeomGrid','scale')],
                       )


    def assemble(self,
                 idx: np.ndarray) -> 'GeomGrid':
        """
        Assemble new grid from index map.

        Parameters
        ----------
        idx : numpy.ndarray of int, shape (:,:,:) or (:,:,:,3)
            GeomGrid of flat indices or coordinate indices.

        Returns
        -------
        updated : damask.GeomGrid
            Updated grid-based geometry.
            Cell count of resulting grid matches shape of index map.
        """
        cells = idx.shape[:3]
        flat = (idx if len(idx.shape)==3 else grid_filters.ravel_index(idx)).flatten(order='F')
        ic = {k: v.reshape((-1,)+v.shape[3:],order='F')[flat]
                  .reshape(cells+v.shape[3:],order='F') for k,v in self.initial_conditions.items()}

        return GeomGrid(material = self.material.flatten(order='F')[flat].reshape(cells,order='F'),
                        size     = self.size,
                        origin   = self.origin,
                        initial_conditions = ic,
                        comments = self.comments+[util.execution_stamp('GeomGrid','assemble')],
                       )


    def renumber(self) -> 'GeomGrid':
        """
        Renumber sorted material indices as 0,...,N-1.

        Returns
        -------
        updated : damask.GeomGrid
            Updated grid-based geometry.
        """
        _,renumbered = np.unique(self.material,return_inverse=True)

        return GeomGrid(material = renumbered.reshape(self.cells),
                        size     = self.size,
                        origin   = self.origin,
                        initial_conditions = self.initial_conditions,
                        comments = self.comments+[util.execution_stamp('GeomGrid','renumber')],
                       )


    def substitute(self,
                   from_material: Union[int,IntSequence],
                   to_material: Union[int,IntSequence]) -> 'GeomGrid':
        """
        Substitute material indices.

        Parameters
        ----------
        from_material : (sequence of) int
            Material indices to be substituted.
        to_material : (sequence of) int
            New material indices.

        Returns
        -------
        updated : damask.GeomGrid
            Updated grid-based geometry.
        """
        material = self.material.copy()
        for f,t in zip(from_material if isinstance(from_material,(Sequence,np.ndarray)) else [from_material],
                       to_material if isinstance(to_material,(Sequence,np.ndarray)) else [to_material]): # ToDo Python 3.10 has strict mode for zip
            material[self.material==f] = t

        return GeomGrid(material = material,
                        size     = self.size,
                        origin   = self.origin,
                        initial_conditions = self.initial_conditions,
                        comments = self.comments+[util.execution_stamp('GeomGrid','substitute')],
                       )


    def sort(self) -> 'GeomGrid':
        """
        Sort material indices such that min(material ID) is located at (0,0,0).

        Returns
        -------
        updated : damask.GeomGrid
            Updated grid-based geometry.
        """
        a = self.material.flatten(order='F')
        from_ma = pd.unique(a)
        sort_idx = np.argsort(from_ma)
        ma = np.unique(a)[sort_idx][np.searchsorted(from_ma,a,sorter = sort_idx)]

        return GeomGrid(material = ma.reshape(self.cells,order='F'),
                        size     = self.size,
                        origin   = self.origin,
                        initial_conditions = self.initial_conditions,
                        comments = self.comments+[util.execution_stamp('GeomGrid','sort')],
                       )


    def clean(self,
              distance: float = np.sqrt(3),
              selection: Optional[IntSequence] = None,
              invert_selection: bool = False,
              periodic: bool = True,
              rng_seed: Optional[NumpyRngSeed] = None) -> 'GeomGrid':
        """
        Smooth grid by selecting most frequent material ID within given stencil at each location.

        Parameters
        ----------
        distance : float, optional
            Voxel distance checked for presence of other materials.
            Defaults to sqrt(3).
        selection : (sequence of) int, optional
            Material IDs to consider. Defaults to all.
        invert_selection : bool, optional
            Consider all material IDs except those in selection. Defaults to False.
        periodic : bool, optional
            Assume grid to be periodic. Defaults to True.
        rng_seed : {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, optional
            A seed to initialize the BitGenerator. Defaults to None.
            If None, then fresh, unpredictable entropy will be pulled from the OS.

        Returns
        -------
        updated : damask.GeomGrid
            Updated grid-based geometry.

        Notes
        -----
        If multiple material IDs are most frequent within a stencil, a random choice is taken.
        """
        def most_frequent(stencil: np.ndarray,
                          selection: Union[None,np.ndarray],
                          rng: np.random.Generator):
            me = stencil[stencil.size//2]
            if selection is None or me in selection:
                unique, counts = np.unique(stencil,return_counts=True)
                return rng.choice(unique[counts==np.max(counts)])
            else:
                return me

        rng = np.random.default_rng(rng_seed)

        d = np.floor(distance).astype(np.int64)
        ext = np.linspace(-d,d,1+2*d,dtype=float),
        xx,yy,zz = np.meshgrid(ext,ext,ext)
        footprint = xx**2+yy**2+zz**2 <= distance**2+distance*1e-8
        selection_ = None if selection is None else \
                     np.setdiff1d(self.material,selection) if invert_selection else \
                     np.intersect1d(self.material,selection)
        material = ndimage.generic_filter(
                                          self.material,
                                          most_frequent,
                                          footprint=footprint,
                                          mode='wrap' if periodic else 'nearest',
                                          extra_keywords=dict(selection=selection_,rng=rng),
                                         ).astype(self.material.dtype)
        return GeomGrid(material = material,
                        size     = self.size,
                        origin   = self.origin,
                        initial_conditions = self.initial_conditions,
                        comments = self.comments+[util.execution_stamp('GeomGrid','clean')],
                       )


    def add_primitive(self,
                      dimension: Union[FloatSequence, IntSequence],
                      center: Union[FloatSequence, IntSequence],
                      exponent: Union[FloatSequence, float],
                      fill: Optional[int] = None,
                      R: Rotation = Rotation(),
                      inverse: bool = False,
                      periodic: bool = True) -> 'GeomGrid':
        """
        Insert a primitive geometric object at a given position.

        The shape of the object is defined as

        .. math::

           |x|^{2^a} + |y|^{2^b} + |z|^{2^c} < 1

        where dimension = (x,y,z) and exponent = (a,b,c).

        Parameters
        ----------
        dimension : sequence of int or float, len (3)
            Dimension (diameter/side length) of the primitive.
            If given as integers, cell centers are addressed.
            If given as floats, physical coordinates are addressed.
        center : sequence of int or float, len (3)
            Center of the primitive.
            If given as integers, cell centers are addressed.
            If given as floats, physical coordinates are addressed.
        exponent : (sequence of) float, len (3)
            Exponents for the three axes.
            0 gives octahedron (ǀxǀ^(2^0) + ǀyǀ^(2^0) + ǀzǀ^(2^0) < 1),
            1 gives sphere     (ǀxǀ^(2^1) + ǀyǀ^(2^1) + ǀzǀ^(2^1) < 1).
        fill : int, optional
            Fill value for primitive. Defaults to material.max()+1.
        R : damask.Rotation, optional
            Rotation of the primitive. Defaults to no rotation.
        inverse : bool, optional
            Retain original materials within primitive and fill outside.
            Defaults to False.
        periodic : bool, optional
            Assume grid to be periodic. Defaults to True.

        Returns
        -------
        updated : damask.GeomGrid
            Updated grid-based geometry.

        Examples
        --------
        Add a sphere at the center.

        >>> import numpy as np
        >>> import damask
        >>> g = damask.GeomGrid(np.zeros([64]*3,int), np.ones(3)*1e-4)
        >>> g.add_primitive(np.ones(3)*5e-5,np.ones(3)*5e-5,1)
        cells:  64 × 64 × 64
        size:   0.0001 × 0.0001 × 0.0001 m³
        origin: 0.0   0.0   0.0 m
        # materials: 2

        Add a cube at the origin.

        >>> import numpy as np
        >>> import damask
        >>> g = damask.GeomGrid(np.zeros([64]*3,int), np.ones(3)*1e-4)
        >>> g.add_primitive(np.ones(3,int)*32,np.zeros(3),np.inf)
        cells:  64 × 64 × 64
        size:   0.0001 × 0.0001 × 0.0001 m³
        origin: 0.0   0.0   0.0 m
        # materials: 2
        """
        # radius and center
        r = np.array(dimension)/2.0*self.size/self.cells if np.issubdtype(np.array(dimension).dtype,np.integer) else \
            np.array(dimension)/2.0
        c = (np.array(center) + .5)*self.size/self.cells if np.issubdtype(np.array(center).dtype,   np.integer) else \
            (np.array(center) - self.origin)

        coords = grid_filters.coordinates0_point(self.cells,self.size,
                                          -(0.5*(self.size + (self.size/self.cells
                                                              if np.issubdtype(np.array(center).dtype,np.integer) else
                                                              0)) if periodic else c))
        coords_rot = R.broadcast_to(tuple(self.cells))@coords

        with np.errstate(all='ignore'):
            mask = np.sum(np.power(np.abs(coords_rot)/r,2.0**np.array(exponent)),axis=-1) > 1.0

        if periodic:                                                                                # translate back to center
            mask = np.roll(mask,((c/self.size-0.5)*self.cells).round().astype(np.int64),(0,1,2))

        return GeomGrid(material = np.where(np.logical_not(mask) if inverse else mask,
                                        self.material,
                                        np.nanmax(self.material)+1 if fill is None else fill),
                        size     = self.size,
                        origin   = self.origin,
                        initial_conditions = self.initial_conditions,
                        comments = self.comments+[util.execution_stamp('GeomGrid','add_primitive')],
                       )


    def vicinity_offset(self,
                        distance: float = np.sqrt(3),
                        offset: Optional[int] = None,
                        selection: Optional[IntSequence] = None,
                        invert_selection: bool = False,
                        periodic: bool = True) -> 'GeomGrid':
        """
        Offset material ID of points in the vicinity of selected (or just other) material IDs.

        Trigger points are variations in material ID, i.e. grain/phase
        boundaries or explicitly given material IDs.

        Parameters
        ----------
        distance : float, optional
            Voxel distance checked for presence of other materials.
            Defaults to sqrt(3).
        offset : int, optional
            Offset (positive or negative) to tag material IDs.
            Defaults to material.max()+1.
        selection : (sequence of) int, optional
            Material IDs that trigger an offset.
            Defaults to any other than own material ID.
        invert_selection : bool, optional
            Consider all material IDs except those in selection.
            Defaults to False.
        periodic : bool, optional
            Assume grid to be periodic. Defaults to True.

        Returns
        -------
        updated : damask.GeomGrid
            Updated grid-based geometry.
        """
        @util.numba_njit_wrapper()
        def tainted_neighborhood(stencil: np.ndarray,
                                 selection: Optional[np.ndarray] = None):
            me = stencil[stencil.size//2]
            if selection is None:
                return np.any(stencil != me)
            elif not len(selection)==0:
                for stencil_item in stencil:
                    for selection_item in selection:
                        if stencil_item==selection_item and selection_item!=me:
                            return True
            return False
        d = np.floor(distance).astype(np.int64)
        ext = np.linspace(-d,d,1+2*d,dtype=float),
        xx,yy,zz = np.meshgrid(ext,ext,ext)
        footprint = xx**2+yy**2+zz**2 <= distance**2+distance*1e-8
        offset_ = np.nanmax(self.material)+1 if offset is None else offset
        selection_ = None if selection is None else \
                     np.setdiff1d(self.material,selection) if invert_selection else \
                     np.intersect1d(self.material,selection)

        mask = ndimage.generic_filter(self.material,
                                      tainted_neighborhood,
                                      footprint=footprint,
                                      mode='wrap' if periodic else 'nearest',
                                      extra_keywords=dict(selection=selection_),
                                     )

        return GeomGrid(material = np.where(mask, self.material + offset_,self.material),
                        size     = self.size,
                        origin   = self.origin,
                        initial_conditions = self.initial_conditions,
                        comments = self.comments+[util.execution_stamp('GeomGrid','vicinity_offset')],
                       )


    def get_grain_boundaries(self,
                             periodic: bool = True,
                             directions: Sequence[str] = 'xyz') -> VTK:
        """
        Create VTK unstructured grid containing grain boundaries.

        Parameters
        ----------
        periodic : bool, optional
            Assume grid to be periodic. Defaults to True.
        directions : (sequence of) {'x', 'y', 'z'}, optional
            Direction(s) along which the boundaries are determined.
            Defaults to 'xyz'.

        Returns
        -------
        grain_boundaries : damask.VTK
            VTK-based geometry of grain boundary network.
        """
        valid = 'xyz'
        if not set(directions).issubset(valid):
            raise ValueError(f'invalid direction "{set(directions).difference(valid)}" specified')

        o = [[0, self.cells[0]+1,           np.prod(self.cells[:2]+1)+self.cells[0]+1, np.prod(self.cells[:2]+1)],
             [0, np.prod(self.cells[:2]+1), np.prod(self.cells[:2]+1)+1,               1],
             [0, 1,                         self.cells[0]+1+1,                         self.cells[0]+1]] # offset for connectivity

        connectivity = []
        for i,d in enumerate(valid):
            if d not in directions: continue
            mask = self.material != np.roll(self.material,1,i)
            for j in [0,1,2]:
                mask = np.concatenate((mask,np.take(mask,[0],j)*(i==j)),j)
            if i == 0 and not periodic: mask[0,:,:] = mask[-1,:,:] = False
            if i == 1 and not periodic: mask[:,0,:] = mask[:,-1,:] = False
            if i == 2 and not periodic: mask[:,:,0] = mask[:,:,-1] = False

            base_nodes = np.argwhere(mask.flatten(order='F')).reshape((-1,1))
            connectivity.append(np.block([base_nodes + o[i][k] for k in range(4)]))

        coords = grid_filters.coordinates0_node(self.cells,self.size,self.origin).reshape((-1,3),order='F')
        return VTK.from_unstructured_grid(coords,np.vstack(connectivity),'QUADRILATERAL')
