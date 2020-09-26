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

    def __init__(self,material,size,origin=[0.0,0.0,0.0],comments=[]):
        """
        New geometry definition from array of material, size, and origin.

        Parameters
        ----------
        material : numpy.ndarray
            Material index array (3D).
        size : list or numpy.ndarray
            Physical size of the geometry in meter.
        origin : list or numpy.ndarray, optional
            Physical origin of the geometry in meter.
        comments : list of str, optional
            Comment lines.

        """
        if len(material.shape) != 3:
            raise ValueError(f'Invalid material shape {material.shape}.')
        elif material.dtype not in np.sctypes['float'] + np.sctypes['int']:
            raise TypeError(f'Invalid material data type {material.dtype}.')
        else:
            self.material = np.copy(material)

            if self.material.dtype in np.sctypes['float'] and \
               np.all(self.material == self.material.astype(int).astype(float)):
                self.material = self.material.astype(int)

        if len(size) != 3 or any(np.array(size) <= 0):
            raise ValueError(f'Invalid size {size}.')
        else:
            self.size = np.array(size)

        if len(origin) != 3:
            raise ValueError(f'Invalid origin {origin}.')
        else:
            self.origin = np.array(origin)

        self.comments = [str(c) for c in comments] if isinstance(comments,list) else [str(comments)]


    def __repr__(self):
        """Basic information on geometry definition."""
        return util.srepr([
               f'grid     a b c:  {util.srepr(self.grid,  " x ")}',
               f'size     x y z:  {util.srepr(self.size,  " x ")}',
               f'origin   x y z:  {util.srepr(self.origin,"   ")}',
               f'# materials:     {self.N_materials}',
               f'max material:    {np.nanmax(self.material)}',
              ])


    def __copy__(self):
        """Copy geometry."""
        return copy.deepcopy(self)


    def copy(self):
        """Copy geometry."""
        return self.__copy__()


    def diff(self,other):
        """
        Report property differences of self relative to other.

        Parameters
        ----------
        other : Geom
            Geometry to compare self against.

        """
        message = []
        if np.any(other.grid != self.grid):
            message.append(util.delete(f'grid     a b c:     {util.srepr(other.grid," x ")}'))
            message.append(util.emph(  f'grid     a b c:     {util.srepr( self.grid," x ")}'))

        if not np.allclose(other.size,self.size):
            message.append(util.delete(f'size     x y z:     {util.srepr(other.size," x ")}'))
            message.append(util.emph(  f'size     x y z:     {util.srepr( self.size," x ")}'))

        if not np.allclose(other.origin,self.origin):
            message.append(util.delete(f'origin   x y z:     {util.srepr(other.origin,"   ")}'))
            message.append(util.emph(  f'origin   x y z:     {util.srepr( self.origin,"   ")}'))

        if other.N_materials != self.N_materials:
            message.append(util.delete(f'# materials:        {other.N_materials}'))
            message.append(util.emph(  f'# materials:        { self.N_materials}'))

        if np.nanmax(other.material) != np.nanmax(self.material):
            message.append(util.delete(f'max material:       {np.nanmax(other.material)}'))
            message.append(util.emph(  f'max material:       {np.nanmax( self.material)}'))

        return util.return_message(message)


    @property
    def grid(self):
        return np.asarray(self.material.shape)


    @property
    def N_materials(self):
        return np.unique(self.material).size


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
            else:
                comments.append(line.strip())

        material = np.empty(grid.prod())                                                      # initialize as flat array
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

        if i != grid.prod():
            raise TypeError(f'Invalid file: expected {grid.prod()} entries, found {i}')

        if not np.any(np.mod(material,1) != 0.0):                                             # no float present
            material = material.astype('int')

        return Geom(material.reshape(grid,order='F'),size,origin,comments)


    @staticmethod
    def load(fname):
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

        return Geom(material = v.get('material').reshape(grid,order='F'),
                    size = bbox[1] - bbox[0],
                    origin = bbox[0],
                    comments=comments)


    @staticmethod
    def _find_closest_seed(seeds, weights, point):
        return np.argmin(np.sum((np.broadcast_to(point,(len(seeds),3))-seeds)**2,axis=1) - weights)

    @staticmethod
    def from_Laguerre_tessellation(grid,size,seeds,weights,material=None,periodic=True):
        """
        Generate geometry from Laguerre tessellation.

        Parameters
        ----------
        grid : int numpy.ndarray of shape (3)
            Number of grid points in x,y,z direction.
        size : list or numpy.ndarray of shape (3)
            Physical size of the geometry in meter.
        seeds : numpy.ndarray of shape (:,3)
            Position of the seed points in meter. All points need to lay within the box.
        weights : numpy.ndarray of shape (seeds.shape[0])
            Weights of the seeds. Setting all weights to 1.0 gives a standard Voronoi tessellation.
        material : numpy.ndarray of shape (seeds.shape[0]), optional
            Material ID of the seeds. Defaults to None, in which case materials are
            consecutively numbered.
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
        material_ = np.array(result.get())

        if periodic:
            material_ = material_.reshape(grid*3)
            material_ = material_[grid[0]:grid[0]*2,grid[1]:grid[1]*2,grid[2]:grid[2]*2]%seeds.shape[0]
        else:
            material_ = material_.reshape(grid)

        return Geom(material = material_+1 if material is None else material[material_],
                    size     = size,
                    comments = util.execution_stamp('Geom','from_Laguerre_tessellation'),
                   )


    @staticmethod
    def from_Voronoi_tessellation(grid,size,seeds,material=None,periodic=True):
        """
        Generate geometry from Voronoi tessellation.

        Parameters
        ----------
        grid : int numpy.ndarray of shape (3)
            Number of grid points in x,y,z direction.
        size : list or numpy.ndarray of shape (3)
            Physical size of the geometry in meter.
        seeds : numpy.ndarray of shape (:,3)
            Position of the seed points in meter. All points need to lay within the box.
        material : numpy.ndarray of shape (seeds.shape[0]), optional
            Material ID of the seeds. Defaults to None, in which case materials are
            consecutively numbered.
        periodic : Boolean, optional
            Perform a periodic tessellation. Defaults to True.

        """
        coords = grid_filters.cell_coord0(grid,size).reshape(-1,3)
        KDTree = spatial.cKDTree(seeds,boxsize=size) if periodic else spatial.cKDTree(seeds)
        devNull,material_ = KDTree.query(coords)

        return Geom(material = (material_+1 if material is None else material[material_]).reshape(grid),
                    size     = size,
                    comments = util.execution_stamp('Geom','from_Voronoi_tessellation'),
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
        header =  [f'{len(self.comments)+4} header'] + self.comments \
                + ['grid   a {} b {} c {}'.format(*self.grid),
                   'size   x {} y {} z {}'.format(*self.size),
                   'origin x {} y {} z {}'.format(*self.origin),
                   'homogenization 1',
                  ]

        grid = self.grid

        if compress is None:
            plain = grid.prod()/self.N_materials < 250
        else:
            plain = not compress

        if plain:
            format_string = '%g' if self.material.dtype in np.sctypes['float'] else \
                            '%{}i'.format(1+int(np.floor(np.log10(np.nanmax(self.material)))))
            np.savetxt(fname,
                       self.material.reshape([grid[0],np.prod(grid[1:])],order='F').T,
                       header='\n'.join(header), fmt=format_string, comments='')
        else:
            try:
                f = open(fname,'w')
            except TypeError:
                f = fname

            compressType = None
            former = start = -1
            reps = 0
            for current in self.material.flatten('F'):
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


    def save(self,fname,compress=True):
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
        v.add(self.material.flatten(order='F'),'material')
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
            Exponents for the three axes.
            0 gives octahedron (ǀxǀ^(2^0) + ǀyǀ^(2^0) + ǀzǀ^(2^0) < 1)
            1 gives sphere     (ǀxǀ^(2^1) + ǀyǀ^(2^1) + ǀzǀ^(2^1) < 1)
        fill : int, optional
            Fill value for primitive. Defaults to material.max() + 1.
        R : damask.Rotation, optional
            Rotation of primitive. Defaults to no rotation.
        inverse : Boolean, optional
            Retain original materials within primitive and fill outside.
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

        fill_ = np.full_like(self.material,np.nanmax(self.material)+1 if fill is None else fill)

        return Geom(material = np.where(np.logical_not(mask) if inverse else mask, self.material,fill_),
                    size     = self.size,
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Geom','add_primitive')],
                   )


    def mirror(self,directions,reflect=False):
        """
        Mirror geometry along given directions.

        Parameters
        ----------
        directions : iterable containing str
            Direction(s) along which the geometry is mirrored.
            Valid entries are 'x', 'y', 'z'.
        reflect : bool, optional
            Reflect (include) outermost layers. Defaults to False.

        """
        valid = ['x','y','z']
        if not set(directions).issubset(valid):
            raise ValueError(f'Invalid direction {set(directions).difference(valid)} specified.')

        limits = [None,None] if reflect else [-2,0]
        mat = self.material.copy()

        if 'x' in directions:
            mat = np.concatenate([mat,mat[limits[0]:limits[1]:-1,:,:]],0)
        if 'y' in directions:
            mat = np.concatenate([mat,mat[:,limits[0]:limits[1]:-1,:]],1)
        if 'z' in directions:
            mat = np.concatenate([mat,mat[:,:,limits[0]:limits[1]:-1]],2)

        return Geom(material = mat,
                    size     = self.size/self.grid*np.asarray(mat.shape),
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Geom','mirror')],
                   )


    def flip(self,directions):
        """
        Flip geometry along given directions.

        Parameters
        ----------
        directions : iterable containing str
            Direction(s) along which the geometry is flipped.
            Valid entries are 'x', 'y', 'z'.

        """
        valid = ['x','y','z']
        if not set(directions).issubset(valid):
            raise ValueError(f'Invalid direction {set(directions).difference(valid)} specified.')

        mat = np.flip(self.material, (valid.index(d) for d in directions if d in valid))

        return Geom(material = mat,
                    size     = self.size,
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Geom','flip')],
                   )


    def scale(self,grid,periodic=True):
        """
        Scale geometry to new grid.

        Parameters
        ----------
        grid : numpy.ndarray of shape (3)
            Number of grid points in x,y,z direction.
        periodic : Boolean, optional
            Assume geometry to be periodic. Defaults to True.

        """
        return Geom(material = ndimage.interpolation.zoom(
                                                           self.material,
                                                           grid/self.grid,
                                                           output=self.material.dtype,
                                                           order=0,
                                                           mode=('wrap' if periodic else 'nearest'),
                                                           prefilter=False
                                                          ),
                    size     = self.size,
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Geom','scale')],
                   )


    def clean(self,stencil=3,selection=None,periodic=True):
        """
        Smooth geometry by selecting most frequent material index within given stencil at each location.

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

        return Geom(material = ndimage.filters.generic_filter(
                                                               self.material,
                                                               mostFrequent,
                                                               size=(stencil if selection is None else stencil//2*2+1,)*3,
                                                               mode=('wrap' if periodic else 'nearest'),
                                                               extra_keywords=dict(selection=selection),
                                                              ).astype(self.material.dtype),
                    size     = self.size,
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Geom','clean')],
                   )


    def renumber(self):
        """Renumber sorted material indices to 1,...,N."""
        renumbered = np.empty(self.grid,dtype=self.material.dtype)
        for i, oldID in enumerate(np.unique(self.material)):
            renumbered = np.where(self.material == oldID, i+1, renumbered)

        return Geom(material = renumbered,
                    size     = self.size,
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Geom','renumber')],
                   )


    def rotate(self,R,fill=None):
        """
        Rotate geometry (pad if required).

        Parameters
        ----------
        R : damask.Rotation
            Rotation to apply to the geometry.
        fill : int or float, optional
            Material index to fill the corners. Defaults to material.max() + 1.

        """
        if fill is None: fill = np.nanmax(self.material) + 1
        dtype = float if np.isnan(fill) or int(fill) != fill or self.material.dtype==np.float else int

        Eulers = R.as_Eulers(degrees=True)
        material_in = self.material.copy()

        # These rotations are always applied in the reference coordinate system, i.e. (z,x,z) not (z,x',z'')
        # see https://www.cs.utexas.edu/~theshark/courses/cs354/lectures/cs354-14.pdf
        for angle,axes in zip(Eulers[::-1], [(0,1),(1,2),(0,1)]):
            material_out = ndimage.rotate(material_in,angle,axes,order=0,
                                                prefilter=False,output=dtype,cval=fill)
            if np.prod(material_in.shape) == np.prod(material_out.shape):
                # avoid scipy interpolation errors for rotations close to multiples of 90°
                material_in = np.rot90(material_in,k=np.rint(angle/90.).astype(int),axes=axes)
            else:
                material_in = material_out

        origin = self.origin-(np.asarray(material_in.shape)-self.grid)*.5 * self.size/self.grid

        return Geom(material = material_in,
                    size     = self.size/self.grid*np.asarray(material_in.shape),
                    origin   = origin,
                    comments = self.comments+[util.execution_stamp('Geom','rotate')],
                   )


    def canvas(self,grid=None,offset=None,fill=None):
        """
        Crop or enlarge/pad geometry.

        Parameters
        ----------
        grid : numpy.ndarray of shape (3)
            Number of grid points in x,y,z direction.
        offset : numpy.ndarray of shape (3)
            Offset (measured in grid points) from old to new geometry [0,0,0].
        fill : int or float, optional
            Material index to fill the background. Defaults to material.max() + 1.

        """
        if offset is None: offset = 0
        if fill is None: fill = np.nanmax(self.material) + 1
        dtype = float if int(fill) != fill or self.material.dtype in np.sctypes['float'] else int

        canvas = np.full(self.grid if grid is None else grid,fill,dtype)

        LL = np.clip( offset,          0,np.minimum(self.grid,     grid+offset))
        UR = np.clip( offset+grid,     0,np.minimum(self.grid,     grid+offset))
        ll = np.clip(-offset,          0,np.minimum(     grid,self.grid-offset))
        ur = np.clip(-offset+self.grid,0,np.minimum(     grid,self.grid-offset))

        canvas[ll[0]:ur[0],ll[1]:ur[1],ll[2]:ur[2]] = self.material[LL[0]:UR[0],LL[1]:UR[1],LL[2]:UR[2]]

        return Geom(material = canvas,
                    size     = self.size/self.grid*np.asarray(canvas.shape),
                    origin   = self.origin+offset*self.size/self.grid,
                    comments = self.comments+[util.execution_stamp('Geom','canvas')],
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
        substituted = self.material.copy()
        for from_ms,to_ms in zip(from_material,to_material):
            substituted[self.material==from_ms] = to_ms

        return Geom(material = substituted,
                    size     = self.size,
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Geom','substitute')],
                   )


    def vicinity_offset(self,vicinity=1,offset=None,trigger=[],periodic=True):
        """
        Offset material index of points in the vicinity of xxx.

        Different from themselves (or listed as triggers) within a given (cubic) vicinity,
        i.e. within the region close to a grain/phase boundary.
        ToDo: use include/exclude as in seeds.from_geom

        Parameters
        ----------
        vicinity : int, optional
            Voxel distance checked for presence of other materials.
            Defaults to 1.
        offset : int, optional
            Offset (positive or negative) to tag material indices,
            defaults to material.max() + 1.
        trigger : list of ints, optional
            List of material indices that trigger a change.
            Defaults to [], meaning that any different neighbor triggers a change.
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

        offset_ = np.nanmax(self.material) if offset is None else offset
        mask = ndimage.filters.generic_filter(self.material,
                                              tainted_neighborhood,
                                              size=1+2*vicinity,
                                              mode='wrap' if periodic else 'nearest',
                                              extra_keywords={'trigger':trigger})

        return Geom(material = np.where(mask, self.material + offset_,self.material),
                    size     = self.size,
                    origin   = self.origin,
                    comments = self.comments+[util.execution_stamp('Geom','vicinity_offset')],
                   )
