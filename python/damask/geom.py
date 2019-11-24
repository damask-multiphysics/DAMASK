import os
from io import StringIO

import numpy as np
from scipy import ndimage
import vtk
from vtk.util import numpy_support

from . import util
from . import version


class Geom():
    """Geometry definition for grid solvers."""

    def __init__(self,microstructure,size,origin=[0.0,0.0,0.0],homogenization=1,comments=[]):
        """
        New geometry definition from array of microstructures and size.

        Parameters
        ----------
        microstructure : numpy.ndarray
            microstructure array (3D)
        size : list or numpy.ndarray
            physical size of the microstructure in meter.
        origin : list or numpy.ndarray, optional
            physical origin of the microstructure in meter.
        homogenization : integer, optional
            homogenization index.
        comments : list of str, optional
            comments lines.

        """
        self.set_microstructure(microstructure)
        self.set_size(size)
        self.set_origin(origin)
        self.set_homogenization(homogenization)
        self.set_comments(comments)
      

    def __repr__(self):
        """Basic information on geometry definition."""
        return util.srepr([
               'grid     a b c:      {}'.format(' x '.join(map(str,self.get_grid  ()))),
               'size     x y z:      {}'.format(' x '.join(map(str,self.get_size  ()))),
               'origin   x y z:      {}'.format('   '.join(map(str,self.get_origin()))),
               'homogenization:      {}'.format(self.get_homogenization()),
               '# microstructures:   {}'.format(len(np.unique(self.microstructure))),
               'max microstructure:  {}'.format(np.nanmax(self.microstructure)),
              ])

    def update(self,microstructure=None,size=None,origin=None,rescale=False):
        """
        Updates microstructure and size.

        Parameters
        ----------
        microstructure : numpy.ndarray, optional
            microstructure array (3D).
        size : list or numpy.ndarray, optional
            physical size of the microstructure in meter.
        origin : list or numpy.ndarray, optional
            physical origin of the microstructure in meter.
        rescale : bool, optional
            ignore size parameter and rescale according to change of grid points.

        """
        grid_old    = self.get_grid()
        size_old    = self.get_size()
        origin_old  = self.get_origin()
        unique_old  = len(np.unique(self.microstructure))
        max_old     = np.nanmax(self.microstructure)
        
        if size is not None and rescale:
            raise ValueError('Either set size explicitly or rescale automatically')

        self.set_microstructure(microstructure)
        self.set_origin(origin)

        if size is not None:
            self.set_size(size)
        elif rescale:
             self.set_size(self.get_grid()/grid_old*self.size)

        message = ['grid     a b c:      {}'.format(' x '.join(map(str,grid_old)))]
        if np.any(grid_old != self.get_grid()):
            message[-1] = util.delete(message[-1])
            message.append(util.emph('grid     a b c:      {}'.format(' x '.join(map(str,self.get_grid())))))

        message.append('size     x y z:      {}'.format(' x '.join(map(str,size_old))))
        if np.any(size_old != self.get_size()):
            message[-1] = util.delete(message[-1])
            message.append(util.emph('size     x y z:      {}'.format(' x '.join(map(str,self.get_size())))))

        message.append('origin   x y z:      {}'.format('   '.join(map(str,origin_old))))
        if np.any(origin_old != self.get_origin()):
            message[-1] = util.delete(message[-1])
            message.append(util.emph('origin   x y z:      {}'.format('   '.join(map(str,self.get_origin())))))

        message.append('homogenization:      {}'.format(self.get_homogenization()))

        message.append('# microstructures:   {}'.format(unique_old))
        if unique_old != len(np.unique(self.microstructure)):
            message[-1] = util.delete(message[-1])
            message.append(util.emph('# microstructures:   {}'.format(len(np.unique(self.microstructure)))))

        message.append('max microstructure:  {}'.format(max_old))
        if max_old != np.nanmax(self.microstructure):
            message[-1] = util.delete(message[-1])
            message.append(util.emph('max microstructure:  {}'.format(np.nanmax(self.microstructure))))
        
        return util.return_message(message)

    def set_comments(self,comments):
        """
        Replaces all existing comments.

        Parameters
        ----------
        comments : list of str
            new comments.

        """
        self.comments = []
        self.add_comments(comments)
      
    def add_comments(self,comments):
        """
        Appends comments to existing comments.

        Parameters
        ----------
        comments : list of str
            new comments.

        """
        self.comments += [str(c) for c in comments] if isinstance(comments,list) else [str(comments)]

    def set_microstructure(self,microstructure):
        """
        Replaces the existing microstructure representation.

        Parameters
        ----------
        microstructure : numpy.ndarray
            microstructure array (3D).

        """
        if microstructure is not None:
            if len(microstructure.shape) != 3:
                raise ValueError('Invalid microstructure shape {}'.format(*microstructure.shape))
            elif microstructure.dtype not in np.sctypes['float'] + np.sctypes['int']:
                raise TypeError('Invalid data type {} for microstructure'.format(microstructure.dtype))
            else:
                self.microstructure = np.copy(microstructure)

    def set_size(self,size):
        """
        Replaces the existing size information.

        Parameters
        ----------
        size : list or numpy.ndarray
            physical size of the microstructure in meter.

        """
        if size is None:
            grid = np.asarray(self.microstructure.shape)
            self.size = grid/np.max(grid)
        else:
            if len(size) != 3 or any(np.array(size)<=0):
                raise ValueError('Invalid size {}'.format(*size))
            else:
                self.size = np.array(size)

    def set_origin(self,origin):
        """
        Replaces the existing origin information.

        Parameters
        ----------
        origin : list or numpy.ndarray
            physical origin of the microstructure in meter

        """
        if origin is not None:
            if len(origin) != 3:
                raise ValueError('Invalid origin {}'.format(*origin))
            else:
                self.origin = np.array(origin)

    def set_homogenization(self,homogenization):
        """
        Replaces the existing homogenization index.

        Parameters
        ----------
        homogenization : integer
            homogenization index

        """
        if homogenization is not None:
            if not isinstance(homogenization,int) or homogenization < 1:
                raise TypeError('Invalid homogenization {}'.format(homogenization))
            else:
                self.homogenization = homogenization


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
        return np.array(self.microstructure.shape)

    def get_homogenization(self):
        """Return the homogenization index."""
        return self.homogenization

    def get_comments(self):
        """Return the comments."""
        return self.comments[:]

    def get_header(self):
        """Return the full header (grid, size, origin, homogenization, comments)."""
        header =  ['{} header'.format(len(self.comments)+4)] + self.comments
        header.append('grid   a {} b {} c {}'.format(*self.get_grid()))
        header.append('size   x {} y {} z {}'.format(*self.get_size()))
        header.append('origin x {} y {} z {}'.format(*self.get_origin()))
        header.append('homogenization {}'.format(self.get_homogenization()))
        return header
        
    @classmethod
    def from_file(cls,fname):
        """
        Reads a geom file.

        Parameters
        ----------
        fname : str or file handle
            geometry file to read.

        """
        with (open(fname) if isinstance(fname,str) else fname) as f:
            f.seek(0)
            header_length,keyword = f.readline().split()[:2]
            header_length = int(header_length)
            content = f.readlines()

        if not keyword.startswith('head') or header_length < 3:
            raise TypeError('Header length information missing or invalid')

        comments = [] 
        for i,line in enumerate(content[:header_length]):
            items = line.lower().strip().split()
            key = items[0] if len(items) > 0 else ''
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
            items = line.split()
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
            raise TypeError('Invalid file: expected {} entries,found {}'.format(grid.prod(),i))
        
        microstructure = microstructure.reshape(grid,order='F')
        if not np.any(np.mod(microstructure.flatten(),1) != 0.0):                                   # no float present
            microstructure = microstructure.astype('int')
        
        return cls(microstructure.reshape(grid),size,origin,homogenization,comments)


    def to_file(self,fname,pack=None):
        """
        Writes a geom file.

        Parameters
        ----------
        fname : str or file handle
            geometry file to write.
        pack : bool, optional
            compress geometry with 'x of y' and 'a to b'.

        """
        header = self.get_header()
        grid   = self.get_grid()
        
        if pack is None:
            plain = grid.prod()/np.unique(self.microstructure).size < 250
        else:
            plain = not pack

        if plain:
            format_string = '%g' if self.microstructure in np.sctypes['float'] else \
                            '%{}i'.format(1+int(np.floor(np.log10(np.nanmax(self.microstructure)))))
            np.savetxt(fname,
                       self.microstructure.reshape([grid[0],np.prod(grid[1:])],order='F').T,
                       header='\n'.join(header), fmt=format_string, comments='')
        else:
            if isinstance(fname,str):
                f = open(fname,'w')
            else:
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
                        f.write('\n'.join(self.get_header())+'\n')
                    elif compressType == '.':
                        f.write('{}\n'.format(former))
                    elif compressType == 'to':
                        f.write('{} to {}\n'.format(start,former))
                    elif compressType == 'of':
                        f.write('{} of {}\n'.format(reps,former))

                    compressType = '.'
                    start = current
                    reps = 1

                former = current

            if compressType == '.':
                f.write('{}\n'.format(former))
            elif compressType == 'to':
                f.write('{} to {}\n'.format(start,former))
            elif compressType == 'of':
                f.write('{} of {}\n'.format(reps,former))

                 
    def to_vtk(self,fname=None):
        """
        Generates vtk file.

        Parameters
        ----------
        fname : str, optional
            vtk file to write. If no file is given, a string is returned.

        """
        grid = self.get_grid() + np.ones(3,dtype=int)
        size = self.get_size()
        origin = self.get_origin()

        coords = [
                  np.linspace(0,size[0],grid[0]) + origin[0],
                  np.linspace(0,size[1],grid[1]) + origin[1],
                  np.linspace(0,size[2],grid[2]) + origin[2]
                 ]
        
        rGrid = vtk.vtkRectilinearGrid()
        coordArray = [vtk.vtkDoubleArray(),vtk.vtkDoubleArray(),vtk.vtkDoubleArray()]

        rGrid.SetDimensions(*grid)
        for d,coord in enumerate(coords):
            for c in coord:
                coordArray[d].InsertNextValue(c)

        rGrid.SetXCoordinates(coordArray[0])
        rGrid.SetYCoordinates(coordArray[1])
        rGrid.SetZCoordinates(coordArray[2])
        
        ms = numpy_support.numpy_to_vtk(num_array=self.microstructure.flatten(order='F'),
                                        array_type=vtk.VTK_INT if self.microstructure.dtype == int else vtk.VTK_FLOAT)
        ms.SetName('microstructure')
        rGrid.GetCellData().AddArray(ms)


        if fname is None:
            writer = vtk.vtkDataSetWriter()
            writer.SetHeader('damask.Geom '+version)
            writer.WriteToOutputStringOn()
        else:
            writer = vtk.vtkXMLRectilinearGridWriter()
            writer.SetCompressorTypeToZLib()
            writer.SetDataModeToBinary()
            
            ext = os.path.splitext(fname)[1]
            if ext == '':
                name = fname + '.' + writer.GetDefaultFileExtension()
            elif ext == writer.GetDefaultFileExtension():
                name = fname
            else:
                raise ValueError("unknown extension {}".format(ext))
            writer.SetFileName(name)
    
        writer.SetInputData(rGrid)
        writer.Write()

        if fname is None: return writer.GetOutputString()

                 
    def show(self):
        """Show raw content (as in file)."""
        f=StringIO()
        self.to_file(f)
        f.seek(0)
        return ''.join(f.readlines())


    def mirror(self,directions,reflect=False):
        """
        Mirror microstructure along given directions.

        Parameters
        ----------
        directions : iterable containing str
            direction(s) along which the microstructure is mirrored. Valid entries are 'x', 'y', 'z'.
        reflect : bool, optional
            reflect (include) outermost layers.

        """
        valid = {'x','y','z'}
        if not all(isinstance(d, str) for d in directions):
            raise TypeError('Directions are not of type str.')
        elif not set(directions).issubset(valid):
            raise ValueError('Invalid direction specified {}'.format(*set(directions).difference(valid)))

        limits = [None,None] if reflect else [-2,0]
        ms = self.get_microstructure()

        if 'z' in directions:
            ms = np.concatenate([ms,ms[:,:,limits[0]:limits[1]:-1]],2)
        if 'y' in directions:
            ms = np.concatenate([ms,ms[:,limits[0]:limits[1]:-1,:]],1)
        if 'x' in directions:
            ms = np.concatenate([ms,ms[limits[0]:limits[1]:-1,:,:]],0)
     
        return self.update(ms,rescale=True)
        #self.add_comments('tbd')


    def clean(self,stencil=3):
        """
        Smooth microstructure by selecting most frequent index within given stencil at each location.

        Parameters
        ----------
        stencil : int, optional
            size of smoothing stencil.

        """
        def mostFrequent(arr):
            unique, inverse = np.unique(arr, return_inverse=True)
            return unique[np.argmax(np.bincount(inverse))]

        return self.update(ndimage.filters.generic_filter(self.microstructure,
                                                          mostFrequent,
                                                          size=(stencil,)*3).astype(self.microstructure.dtype))
        #self.add_comments('tbd')
