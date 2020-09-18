import copy
import multiprocessing as mp
from functools import partial

import numpy as np
from scipy import ndimage,spatial

from . import environment
from . import Rotation
from . import VTK
from . import util
from . import grid_filters


class Geom:
    """Geometry definition for grid solvers."""

    def __init__(self,microstructure,size,origin=[0.0,0.0,0.0],homogenization=1,comments=[]):
        """
        New geometry definition from array of microstructures and size.

        Parameters
        ----------
        microstructure : numpy.ndarray
            Microstructure array (3D)
        size : list or numpy.ndarray
            Physical size of the microstructure in meter.
        origin : list or numpy.ndarray, optional
            Physical origin of the microstructure in meter.
        homogenization : int, optional
            Homogenization index.
        comments : list of str, optional
            Comment lines.

        """
        self.set_microstructure(microstructure)
        self.set_size(size)
        self.set_origin(origin)
        self.set_homogenization(homogenization)
        self.set_comments(comments)


    def __repr__(self):
        """Basic information on geometry definition."""
        return util.srepr([
               f'grid     a b c:     {util.srepr(self.get_grid  ()," x ")}',
               f'size     x y z:     {util.srepr(self.get_size  ()," x ")}',
               f'origin   x y z:     {util.srepr(self.get_origin(),"   ")}',
               f'# materialpoints:   {self.N_microstructure}',
               f'max materialpoint:  {np.nanmax(self.microstructure)}',
              ])


    def __copy__(self):
        """Copy geometry."""
        return copy.deepcopy(self)


    def copy(self):
        """Copy geometry."""
        return self.__copy__()


    def duplicate(self,microstructure=None,size=None,origin=None,comments=None,autosize=False):
        """
        Create a duplicate having updated microstructure, size, and origin.

        Parameters
        ----------
        microstructure : numpy.ndarray, optional
            Microstructure array (3D).
        size : list or numpy.ndarray, optional
            Physical size of the microstructure in meter.
        origin : list or numpy.ndarray, optional
            Physical origin of the microstructure in meter.
        comments : list of str, optional
            Comment lines.
        autosize : bool, optional
            Ignore size parameter and rescale according to change of grid points.

        """
        if size is not None and autosize:
            raise ValueError('Auto-sizing conflicts with explicit size parameter.')

        grid_old = self.get_grid()
        dup = self.copy()
        dup.set_microstructure(microstructure)
        dup.set_origin(origin)

        if comments is not None:
            dup.set_comments(comments)

        if size is not None:
            dup.set_size(size)
        elif autosize:
            dup.set_size(dup.get_grid()/grid_old*self.get_size())

        return dup


    def diff(self,other):
        """
        Report property differences of self relative to other.

        Parameters
        ----------
        other : Geom
            Geometry to compare self against.

        """
        message = []
        if np.any(other.get_grid() != self.get_grid()):
            message.append(util.delete(f'grid     a b c:     {util.srepr(other.get_grid()," x ")}'))
            message.append(util.emph(  f'grid     a b c:     {util.srepr( self.get_grid()," x ")}'))

        if np.any(other.get_size() != self.get_size()):
            message.append(util.delete(f'size     x y z:     {util.srepr(other.get_size()," x ")}'))
            message.append(util.emph(  f'size     x y z:     {util.srepr( self.get_size()," x ")}'))

        if np.any(other.get_origin() != self.get_origin()):
            message.append(util.delete(f'origin   x y z:     {util.srepr(other.get_origin(),"   ")}'))
            message.append(util.emph(  f'origin   x y z:     {util.srepr( self.get_origin(),"   ")}'))

        if other.N_microstructure != self.N_microstructure:
            message.append(util.delete(f'# materialpoints:   {other.N_microstructure}'))
            message.append(util.emph(  f'# materialpoints:   { self.N_microstructure}'))

        if np.nanmax(other.microstructure) != np.nanmax(self.microstructure):
            message.append(util.delete(f'max materialpoint:  {np.nanmax(other.microstructure)}'))
            message.append(util.emph(  f'max materialpoint:  {np.nanmax( self.microstructure)}'))

        return util.return_message(message)


    def set_comments(self,comments):
        """
        Replace all existing comments.

        Parameters
        ----------
        comments : list of str
            All comments.

        """
        self.comments = []
        self.add_comments(comments)


    def add_comments(self,comments):
        """
        Append comments to existing comments.

        Parameters
        ----------
        comments : list of str
            New comments.

        """
        self.comments += [str(c) for c in comments] if isinstance(comments,list) else [str(comments)]


    def set_microstructure(self,microstructure):
        """
        Replace the existing microstructure representation.

        The complete microstructure is replaced (indcluding grid definition),
        unless a masked array is provided in which case the grid dimensions
        need to match and masked entries are not replaced.

        Parameters
        ----------
        microstructure : numpy.ndarray or numpy.ma.core.MaskedArray of shape (:,:,:)
            Microstructure indices.

        """
        if microstructure is not None:
            if isinstance(microstructure,np.ma.core.MaskedArray):
                self.microstructure = np.where(microstructure.mask,
                                               self.microstructure,microstructure.data)
            else:
                self.microstructure = np.copy(microstructure)

        if         self.microstructure.dtype in np.sctypes['float'] and \
            np.all(self.microstructure == self.microstructure.astype(int).astype(float)):
            self.microstructure = self.microstructure.astype(int)

        if len(self.microstructure.shape) != 3:
            raise ValueError(f'Invalid microstructure shape {microstructure.shape}')
        elif self.microstructure.dtype not in np.sctypes['float'] + np.sctypes['int']:
            raise TypeError(f'Invalid microstructure data type {microstructure.dtype}')


    def set_size(self,size):
        """
        Replace the existing size information.

        Parameters
        ----------
        size : list or numpy.ndarray
            Physical size of the microstructure in meter.

        """
        if size is not None:
            if len(size) != 3 or any(np.array(size) <= 0):
                raise ValueError(f'Invalid size {size}')
            else:
                self.size = np.array(size)


    def set_origin(self,origin):
        """
        Replace the existing origin information.

        Parameters
        ----------
        origin : list or numpy.ndarray
            Physical origin of the microstructure in meter.

        """
        if origin is not None:
            if len(origin) != 3:
                raise ValueError(f'Invalid origin {origin}')
            else:
                self.origin = np.array(origin)


    def set_homogenization(self,homogenization):
        """
        Replace the existing homogenization index.

        Parameters
        ----------
        homogenization : int
            Homogenization index.

        """
        if homogenization is not None:
            if not isinstance(homogenization,int) or homogenization < 1:
                raise TypeError(f'Invalid homogenization {homogenization}.')
            else:
                self.homogenization = homogenization


    @property
    def grid(self):
        return self.get_grid()


    @property
    def N_microstructure(self):
        return np.unique(self.microstructure).size


    def get_microstructure(self):
        """Return the microstructure representation."""
        return np.copy(self.microstructure)


    def get_size(self):
        """Return the physical size in meter."""
        return np.copy(self.size)


    def get_origin(self):
        """Return the origin in meter."""
        return np.copy(self.origin)


    def get_grid(self):
        """Return the grid discretization."""
        return np.asarray(self.microstructure.shape)


    def get_homogenization(self):
        """Return the homogenization index."""
        return self.homogenization


    def get_comments(self):
        """Return the comments."""
        return self.comments[:]


    @staticmethod
    def load_ASCII(fname):
        """
        Read a geom file.

        Parameters
        ----------
        fname : str or file handle
            Geometry file to read.

        """
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
            raise TypeError('Header length information missing or invalid')

        content = f.readlines()

        comments = []
        for i,line in enumerate(content[:header_length]):
            items = line.split('#')[0].lower().strip().split()
            key = items[0] if items else ''
            if   key == 'grid':
                grid   = np.array([  int(dict(zip(items[1::2],items[2::2]))[i]) for i in ['a','b','c']])
            elif key == 'size':
                size   = np.array([float(dict(zip(items[1::2],items[2::2]))[i]) for i in ['x','y','z']])
            elif key == 'origin':
                origin = np.array([float(dict(zip(items[1::2],items[2::2]))[i]) for i in ['x','y','z']])
            elif key == 'homogenization':
                homogenization = int(items[1])
            else:
                comments.append(line.strip())

        microstructure = np.empty(grid.prod())                                                      # initialize as flat array
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
            microstructure[i:i+len(items)] = items
            i += len(items)

        if i != grid.prod():
            raise TypeError(f'Invalid file: expected {grid.prod()} entries, found {i}')

        if not np.any(np.mod(microstructure,1) != 0.0):                                             # no float present
            microstructure = microstructure.astype('int')

        return Geom(microstructure.reshape(grid,order='F'),size,origin,homogenization,comments)


    @staticmethod
    def load_vtr(fname):
        """
        Read a VTK rectilinear grid.

        Parameters
        ----------
        fname : str or or pathlib.Path
            Geometry file to read.
            Valid extension is .vtr, it will be appended if not given.

        """
        v = VTK.load(fname if str(fname).endswith('.vtr') else str(fname)+'.vtr')
        comments = v.get_comments()
        grid = np.array(v.vtk_data.GetDimensions())-1
        bbox = np.array(v.vtk_data.GetBounds()).reshape(3,2).T
        size = bbox[1] - bbox[0]

        return Geom(v.get('materialpoint').reshape(grid,order='F'),size,bbox[0],comments=comments)


    @staticmethod
    def _find_closest_seed(seeds, weights, point):
        return np.argmin(np.sum((np.broadcast_to(point,(len(seeds),3))-seeds)**2,axis=1) - weights)

    @staticmethod
    def from_Laguerre_tessellation(grid,size,seeds,weights,periodic=True):
        """
        Generate geometry from Laguerre tessellation.

        Parameters
        ----------
        grid : int numpy.ndarray of shape (3)
            Number of grid points in x,y,z direction.
        size : list or numpy.ndarray of shape (3)
            Physical size of the microstructure in meter.
        seeds : numpy.ndarray of shape (:,3)
            Position of the seed points in meter. All points need to lay within the box.
        weights : numpy.ndarray of shape (seeds.shape[0])
            Weights of the seeds. Setting all weights to 1.0 gives a standard Voronoi tessellation.
        periodic : Boolean, optional
            Perform a periodic tessellation. Defaults to True.

        """
        if periodic:
            weights_p = np.tile(weights,27)                                                         # Laguerre weights (1,2,3,1,2,3,...,1,2,3)
            seeds_p = np.vstack((seeds  -np.array([size[0],0.,0.]),seeds,  seeds  +np.array([size[0],0.,0.])))
            seeds_p = np.vstack((seeds_p-np.array([0.,size[1],0.]),seeds_p,seeds_p+np.array([0.,size[1],0.])))
            seeds_p = np.vstack((seeds_p-np.array([0.,0.,size[2]]),seeds_p,seeds_p+np.array([0.,0.,size[2]])))
            coords  = grid_filters.cell_coord0(grid*3,size*3,-size).reshape(-1,3)
        else:
            weights_p = weights
            seeds_p   = seeds
            coords    = grid_filters.cell_coord0(grid,size).reshape(-1,3)

        pool = mp.Pool(processes = int(environment.options['DAMASK_NUM_THREADS']))
        result = pool.map_async(partial(Geom._find_closest_seed,seeds_p,weights_p), [coord for coord in coords])
        pool.close()
        pool.join()
        microstructure = np.array(result.get())

        if periodic:
            microstructure = microstructure.reshape(grid*3)
            microstructure = microstructure[grid[0]:grid[0]*2,grid[1]:grid[1]*2,grid[2]:grid[2]*2]%seeds.shape[0]
        else:
            microstructure = microstructure.reshape(grid)

        return Geom(microstructure+1,size,homogenization=1,
                    comments=util.execution_stamp('Geom','from_Laguerre_tessellation'),
                   )


    @staticmethod
    def from_Voronoi_tessellation(grid,size,seeds,periodic=True):
        """
        Generate geometry from Voronoi tessellation.

        Parameters
        ----------
        grid : int numpy.ndarray of shape (3)
            Number of grid points in x,y,z direction.
        size : list or numpy.ndarray of shape (3)
            Physical size of the microstructure in meter.
        seeds : numpy.ndarray of shape (:,3)
            Position of the seed points in meter. All points need to lay within the box.
        periodic : Boolean, optional
            Perform a periodic tessellation. Defaults to True.

        """
        coords = grid_filters.cell_coord0(grid,size).reshape(-1,3)
        KDTree = spatial.cKDTree(seeds,boxsize=size) if periodic else spatial.cKDTree(seeds)
        devNull,microstructure = KDTree.query(coords)

        return Geom(microstructure.reshape(grid)+1,size,homogenization=1,
                    comments=util.execution_stamp('Geom','from_Voronoi_tessellation'),
                   )


    def save_ASCII(self,fname,compress=None):
        """
        Writes a geom file.

        Parameters
        ----------
        fname : str or file handle
            Geometry file to write with extension '.geom'.
        compress : bool, optional
            Compress geometry with 'x of y' and 'a to b'.

        """
        header =  [f'{len(self.comments)+4} header'] + self.comments
        header.append('grid   a {} b {} c {}'.format(*self.get_grid()))
        header.append('size   x {} y {} z {}'.format(*self.get_size()))
        header.append('origin x {} y {} z {}'.format(*self.get_origin()))
        header.append(f'homogenization {self.get_homogenization()}')

        grid = self.get_grid()

        if compress is None:
            plain = grid.prod()/self.N_microstructure < 250
        else:
            plain = not compress

        if plain:
            format_string = '%g' if self.microstructure.dtype in np.sctypes['float'] else \
                            '%{}i'.format(1+int(np.floor(np.log10(np.nanmax(self.microstructure)))))
            np.savetxt(fname,
                       self.microstructure.reshape([grid[0],np.prod(grid[1:])],order='F').T,
                       header='\n'.join(header), fmt=format_string, comments='')
        else:
            try:
                f = open(fname,'w')
            except TypeError:
                f = fname

            compressType = None
            former = start = -1
            reps = 0
            for current in self.microstructure.flatten('F'):
                if abs(current - former) == 1 and (start - current) == reps*(former - current):
                    compressType = 'to'
                    reps += 1
                elif current == former and start == former:
                    compressType = 'of'
                    reps += 1
                else:
                    if   compressType is None:
                        f.write('\n'.join(header)+'\n')
                    elif compressType == '.':
                        f.write(f'{former}\n')
                    elif compressType == 'to':
                        f.write(f'{start} to {former}\n')
                    elif compressType == 'of':
                        f.write(f'{reps} of {former}\n')

                    compressType = '.'
                    start = current
                    reps = 1

                former = current

            if compressType == '.':
                f.write(f'{former}\n')
            elif compressType == 'to':
                f.write(f'{start} to {former}\n')
            elif compressType == 'of':
                f.write(f'{reps} of {former}\n')


    def save_vtr(self,fname,compress=True):
        """
        Generates vtk rectilinear grid.

        Parameters
        ----------
        fname : str, optional
            Filename to write. If no file is given, a string is returned.
            Valid extension is .vtr, it will be appended if not given.
        compress : bool, optional
            Compress with zlib algorithm. Defaults to True.

        """
        v = VTK.from_rectilinearGrid(self.grid,self.size,self.origin)
        v.add(self.microstructure.flatten(order='F'),'materialpoint')
        v.add_comments(self.comments)

        v.save(fname if str(fname).endswith('.vtr') else str(fname)+'.vtr',parallel=False,compress=compress)


    def show(self):
        """Show on screen."""
        v = VTK.from_rectilinearGrid(self.grid,self.size,self.origin)
        v.show()


    def add_primitive(self,dimension,center,exponent,
                      fill=None,R=Rotation(),inverse=False,periodic=True):
        """
        Inserts a primitive geometric object at a given position.

        Parameters
        ----------
        dimension : int or float numpy.ndarray of shape(3)
            Dimension (diameter/side length) of the primitive. If given as
            integers, grid point locations (cell centers) are addressed.
            If given as floats, coordinates are addressed.
        center : int or float numpy.ndarray of shape(3)
            Center of the primitive. If given as integers, grid point
            locations (cell centers) are addressed.
            If given as floats, coordinates are addressed.
        exponent : numpy.ndarray of shape(3) or float
            Exponents for the three axis.
            0 gives octahedron (|x|^(2^0) + |y|^(2^0) + |z|^(2^0) < 1)
            1 gives a sphere (|x|^(2^1) + |y|^(2^1) + |z|^(2^1) < 1)
        fill : int, optional
            Fill value for primitive. Defaults to microstructure.max() + 1.
        R : damask.Rotation, optional
            Rotation of primitive. Defaults to no rotation.
        inverse : Boolean, optional
            Retain original microstructure within primitive and fill outside.
            Defaults to False.
        periodic : Boolean, optional
            Repeat primitive over boundaries. Defaults to True.

        """
        # normalized 'radius' and center
        r = np.array(dimension)/self.grid/2.0 if np.array(dimension).dtype in np.sctypes['int'] else \
            np.array(dimension)/self.size/2.0
        c = (np.array(center) + .5)/self.grid if np.array(center).dtype in np.sctypes['int'] else \
            (np.array(center) - self.origin)/self.size

        coords = grid_filters.cell_coord0(self.grid,np.ones(3)) \
               - ((np.ones(3)-(1./self.grid if np.array(center).dtype in np.sctypes['int'] else 0))*0.5 if periodic else c)  # periodic center is always at CoG
        coords_rot = R.broadcast_to(tuple(self.grid))@coords

        with np.errstate(all='ignore'):
            mask = np.where(np.sum(np.power(coords_rot/r,2.0**exponent),axis=-1) > 1.0,True,False)

        if periodic:                                                                                # translate back to center
            mask = np.roll(mask,((c-np.ones(3)*.5)*self.grid).astype(int),(0,1,2))

        fill_ = np.full_like(self.microstructure,np.nanmax(self.microstructure)+1 if fill is None else fill)
        ms = np.ma.MaskedArray(fill_,np.logical_not(mask) if inverse else mask)

        return self.duplicate(ms,
                              comments=self.get_comments()+[util.execution_stamp('Geom','add_primitive')],
                             )


    def mirror(self,directions,reflect=False):
        """
        Mirror microstructure along given directions.

        Parameters
        ----------
        directions : iterable containing str
            Direction(s) along which the microstructure is mirrored.
            Valid entries are 'x', 'y', 'z'.
        reflect : bool, optional
            Reflect (include) outermost layers. Defaults to False.

        """
        valid = ['x','y','z']
        if not set(directions).issubset(valid):
            raise ValueError(f'Invalid direction {set(directions).difference(valid)} specified.')

        limits = [None,None] if reflect else [-2,0]
        ms = self.get_microstructure()

        if 'z' in directions:
            ms = np.concatenate([ms,ms[:,:,limits[0]:limits[1]:-1]],2)
        if 'y' in directions:
            ms = np.concatenate([ms,ms[:,limits[0]:limits[1]:-1,:]],1)
        if 'x' in directions:
            ms = np.concatenate([ms,ms[limits[0]:limits[1]:-1,:,:]],0)

        return self.duplicate(ms,
                              comments=self.get_comments()+[util.execution_stamp('Geom','mirror')],
                              autosize=True)


    def flip(self,directions):
        """
        Flip microstructure along given directions.

        Parameters
        ----------
        directions : iterable containing str
            Direction(s) along which the microstructure is flipped.
            Valid entries are 'x', 'y', 'z'.

        """
        valid = ['x','y','z']
        if not set(directions).issubset(valid):
            raise ValueError(f'Invalid direction {set(directions).difference(valid)} specified.')

        ms = np.flip(self.microstructure, (valid.index(d) for d in directions if d in valid))

        return self.duplicate(ms,
                              comments=self.get_comments()+[util.execution_stamp('Geom','flip')],
                             )


    def scale(self,grid,periodic=True):
        """
        Scale microstructure to new grid.

        Parameters
        ----------
        grid : numpy.ndarray of shape (3)
            Number of grid points in x,y,z direction.
        periodic : Boolean, optional
            Assume geometry to be periodic. Defaults to True.

        """
        return self.duplicate(ndimage.interpolation.zoom(
                                                         self.microstructure,
                                                         grid/self.get_grid(),
                                                         output=self.microstructure.dtype,
                                                         order=0,
                                                         mode=('wrap' if periodic else 'nearest'),
                                                         prefilter=False
                                                        ),
                              comments=self.get_comments()+[util.execution_stamp('Geom','scale')],
                             )


    def clean(self,stencil=3,selection=None,periodic=True):
        """
        Smooth microstructure by selecting most frequent index within given stencil at each location.

        Parameters
        ----------
        stencil : int, optional
            Size of smoothing stencil.
        selection : list, optional
            Field values that can be altered. Defaults to all.
        periodic : Boolean, optional
            Assume geometry to be periodic. Defaults to True.

        """
        def mostFrequent(arr,selection=None):
            me = arr[arr.size//2]
            if selection is None or me in selection:
                unique, inverse = np.unique(arr, return_inverse=True)
                return unique[np.argmax(np.bincount(inverse))]
            else:
                return me

        return self.duplicate(ndimage.filters.generic_filter(
                                                             self.microstructure,
                                                             mostFrequent,
                                                             size=(stencil if selection is None else stencil//2*2+1,)*3,
                                                             mode=('wrap' if periodic else 'nearest'),
                                                             extra_keywords=dict(selection=selection),
                                                            ).astype(self.microstructure.dtype),
                              comments=self.get_comments()+[util.execution_stamp('Geom','clean')],
                             )


    def renumber(self):
        """Renumber sorted microstructure indices to 1,...,N."""
        renumbered = np.empty(self.get_grid(),dtype=self.microstructure.dtype)
        for i, oldID in enumerate(np.unique(self.microstructure)):
            renumbered = np.where(self.microstructure == oldID, i+1, renumbered)

        return self.duplicate(renumbered,
                              comments=self.get_comments()+[util.execution_stamp('Geom','renumber')],
                             )


    def rotate(self,R,fill=None):
        """
        Rotate microstructure (pad if required).

        Parameters
        ----------
        R : damask.Rotation
            Rotation to apply to the microstructure.
        fill : int or float, optional
            Microstructure index to fill the corners. Defaults to microstructure.max() + 1.

        """
        if fill is None: fill = np.nanmax(self.microstructure) + 1
        dtype = float if np.isnan(fill) or int(fill) != fill or self.microstructure.dtype==np.float else int

        Eulers = R.as_Eulers(degrees=True)
        microstructure_in = self.get_microstructure()

        # These rotations are always applied in the reference coordinate system, i.e. (z,x,z) not (z,x',z'')
        # see https://www.cs.utexas.edu/~theshark/courses/cs354/lectures/cs354-14.pdf
        for angle,axes in zip(Eulers[::-1], [(0,1),(1,2),(0,1)]):
            microstructure_out = ndimage.rotate(microstructure_in,angle,axes,order=0,
                                                prefilter=False,output=dtype,cval=fill)
            if np.prod(microstructure_in.shape) == np.prod(microstructure_out.shape):
                # avoid scipy interpolation errors for rotations close to multiples of 90°
                microstructure_in = np.rot90(microstructure_in,k=np.rint(angle/90.).astype(int),axes=axes)
            else:
                microstructure_in = microstructure_out

        origin = self.origin-(np.asarray(microstructure_in.shape)-self.grid)*.5 * self.size/self.grid

        return self.duplicate(microstructure_in,
                              origin=origin,
                              comments=self.get_comments()+[util.execution_stamp('Geom','rotate')],
                              autosize=True,
                             )


    def canvas(self,grid=None,offset=None,fill=None):
        """
        Crop or enlarge/pad microstructure.

        Parameters
        ----------
        grid : numpy.ndarray of shape (3)
            Number of grid points in x,y,z direction.
        offset : numpy.ndarray of shape (3)
            Offset (measured in grid points) from old to new microstructure[0,0,0].
        fill : int or float, optional
            Microstructure index to fill the background. Defaults to microstructure.max() + 1.

        """
        if offset is None: offset = 0
        if fill is None: fill = np.nanmax(self.microstructure) + 1
        dtype = float if int(fill) != fill or self.microstructure.dtype in np.sctypes['float'] else int

        canvas = np.full(self.grid if grid is None else grid,fill,dtype)

        LL = np.clip( offset,          0,np.minimum(self.grid,     grid+offset))
        UR = np.clip( offset+grid,     0,np.minimum(self.grid,     grid+offset))
        ll = np.clip(-offset,          0,np.minimum(     grid,self.grid-offset))
        ur = np.clip(-offset+self.grid,0,np.minimum(     grid,self.grid-offset))

        canvas[ll[0]:ur[0],ll[1]:ur[1],ll[2]:ur[2]] = self.microstructure[LL[0]:UR[0],LL[1]:UR[1],LL[2]:UR[2]]

        return self.duplicate(canvas,
                              origin=self.origin+offset*self.size/self.grid,
                              comments=self.get_comments()+[util.execution_stamp('Geom','canvas')],
                              autosize=True,
                             )


    def substitute(self,from_microstructure,to_microstructure):
        """
        Substitute microstructure indices.

        Parameters
        ----------
        from_microstructure : iterable of ints
            Microstructure indices to be substituted.
        to_microstructure : iterable of ints
            New microstructure indices.

        """
        substituted = self.get_microstructure()
        for from_ms,to_ms in zip(from_microstructure,to_microstructure):
            substituted[self.microstructure==from_ms] = to_ms

        return self.duplicate(substituted,
                              comments=self.get_comments()+[util.execution_stamp('Geom','substitute')],
                             )


    def vicinity_offset(self,vicinity=1,offset=None,trigger=[],periodic=True):
        """
        Offset microstructure index of points in the vicinity of xxx.

        Different from themselves (or listed as triggers) within a given (cubic) vicinity,
        i.e. within the region close to a grain/phase boundary.
        ToDo: use include/exclude as in seeds.from_geom

        Parameters
        ----------
        vicinity : int, optional
            Voxel distance checked for presence of other microstructure.
            Defaults to 1.
        offset : int, optional
            Offset (positive or negative) to tag microstructure indices,
            defaults to microstructure.max() + 1.
        trigger : list of ints, optional
            List of microstructure indices triggering a change.
            Defaults to [], meaning that different neigboors trigger a change.
        periodic : Boolean, optional
            Assume geometry to be periodic. Defaults to True.

        """
        def tainted_neighborhood(stencil,trigger):

            me = stencil[stencil.shape[0]//2]
            if len(trigger) == 0:
                return np.any(stencil != me)
            if me in trigger:
                trigger = set(trigger)
                trigger.remove(me)
                trigger = list(trigger)
            return np.any(np.in1d(stencil,np.array(trigger)))

        offset_ = np.nanmax(self.microstructure) if offset is None else offset
        mask = ndimage.filters.generic_filter(self.microstructure,
                                              tainted_neighborhood,
                                              size=1+2*vicinity,
                                              mode='wrap' if periodic else 'nearest',
                                              extra_keywords={'trigger':trigger})
        microstructure = np.ma.MaskedArray(self.microstructure + offset_, np.logical_not(mask))

        return self.duplicate(microstructure,
                              comments=self.get_comments()+[util.execution_stamp('Geom','vicinity_offset')],
                             )
