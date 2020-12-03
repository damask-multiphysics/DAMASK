"""
Filters for operations on regular grids.

Notes
-----
The grids are defined as (x,y,z,...) where x is fastest and z is slowest.
This convention is consistent with the geom file format.
When converting to/from a plain list (e.g. storage in ASCII table),
the following operations are required for tensorial data:

D3 = D1.reshape(cells+(-1,),order='F').reshape(cells+(3,3))
D1 = D3.reshape(cells+(-1,)).reshape(-1,9,order='F')

"""
from scipy import spatial as _spatial
import numpy as _np

def _ks(size,cells,first_order=False):
    """
    Get wave numbers operator.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        physical size of the periodic field.
    cells : numpy.ndarray of shape (3)
        number of cells.
    first_order : bool, optional
        correction for first order derivatives, defaults to False.

    """
    k_sk = _np.where(_np.arange(cells[0])>cells[0]//2,_np.arange(cells[0])-cells[0],_np.arange(cells[0]))/size[0]
    if cells[0]%2 == 0 and first_order: k_sk[cells[0]//2] = 0                                       # Nyquist freq=0 for even cells (Johnson, MIT, 2011)

    k_sj = _np.where(_np.arange(cells[1])>cells[1]//2,_np.arange(cells[1])-cells[1],_np.arange(cells[1]))/size[1]
    if cells[1]%2 == 0 and first_order: k_sj[cells[1]//2] = 0                                       # Nyquist freq=0 for even cells (Johnson, MIT, 2011)

    k_si = _np.arange(cells[2]//2+1)/size[2]

    return _np.stack(_np.meshgrid(k_sk,k_sj,k_si,indexing = 'ij'), axis=-1)


def curl(size,field):
    """
    Calculate curl of a vector or tensor field in Fourier space.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        physical size of the periodic field.
    field : numpy.ndarray of shape (:,:,:,3) or (:,:,:,3,3)
        periodic field of which the curl is calculated.

    """
    n = _np.prod(field.shape[3:])
    k_s = _ks(size,field.shape[:3],True)

    e = _np.zeros((3, 3, 3))
    e[0, 1, 2] = e[1, 2, 0] = e[2, 0, 1] = +1.0                                                     # Levi-Civita symbol
    e[0, 2, 1] = e[2, 1, 0] = e[1, 0, 2] = -1.0

    field_fourier = _np.fft.rfftn(field,axes=(0,1,2))
    curl_ = (_np.einsum('slm,ijkl,ijkm ->ijks', e,k_s,field_fourier)*2.0j*_np.pi if n == 3 else     # vector, 3   -> 3
             _np.einsum('slm,ijkl,ijknm->ijksn',e,k_s,field_fourier)*2.0j*_np.pi)                   # tensor, 3x3 -> 3x3

    return _np.fft.irfftn(curl_,axes=(0,1,2),s=field.shape[:3])


def divergence(size,field):
    """
    Calculate divergence of a vector or tensor field in Fourier space.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        physical size of the periodic field.
    field : numpy.ndarray of shape (:,:,:,3) or (:,:,:,3,3)
        periodic field of which the divergence is calculated.

    """
    n = _np.prod(field.shape[3:])
    k_s = _ks(size,field.shape[:3],True)

    field_fourier = _np.fft.rfftn(field,axes=(0,1,2))
    div_ = (_np.einsum('ijkl,ijkl ->ijk', k_s,field_fourier)*2.0j*_np.pi if n == 3 else             # vector, 3   -> 1
            _np.einsum('ijkm,ijklm->ijkl',k_s,field_fourier)*2.0j*_np.pi)                           # tensor, 3x3 -> 3

    return _np.fft.irfftn(div_,axes=(0,1,2),s=field.shape[:3])


def gradient(size,field):
    """
    Calculate gradient of a scalar or vector field in Fourier space.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        physical size of the periodic field.
    field : numpy.ndarray of shape (:,:,:,1) or (:,:,:,3)
        periodic field of which the gradient is calculated.

    """
    n = _np.prod(field.shape[3:])
    k_s = _ks(size,field.shape[:3],True)

    field_fourier = _np.fft.rfftn(field,axes=(0,1,2))
    grad_ = (_np.einsum('ijkl,ijkm->ijkm', field_fourier,k_s)*2.0j*_np.pi if n == 1 else            # scalar, 1 -> 3
             _np.einsum('ijkl,ijkm->ijklm',field_fourier,k_s)*2.0j*_np.pi)                          # vector, 3 -> 3x3

    return _np.fft.irfftn(grad_,axes=(0,1,2),s=field.shape[:3])


def cell_coord0(cells,size,origin=_np.zeros(3)):
    """
    Cell center positions (undeformed).

    Parameters
    ----------
    cells : numpy.ndarray of shape (3)
        number of cells.
    size : numpy.ndarray of shape (3)
        physical size of the periodic field.
    origin : numpy.ndarray, optional
        physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    """
    start = origin        + size/cells*.5
    end   = origin + size - size/cells*.5

    return _np.stack(_np.meshgrid(_np.linspace(start[0],end[0],cells[0]),
                                  _np.linspace(start[1],end[1],cells[1]),
                                  _np.linspace(start[2],end[2],cells[2]),indexing = 'ij'),
                     axis = -1)


def cell_displacement_fluct(size,F):
    """
    Cell center displacement field from fluctuation part of the deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.

    """
    integrator = 0.5j*size/_np.pi

    k_s = _ks(size,F.shape[:3],False)
    k_s_squared = _np.einsum('...l,...l',k_s,k_s)
    k_s_squared[0,0,0] = 1.0

    displacement = -_np.einsum('ijkml,ijkl,l->ijkm',
                              _np.fft.rfftn(F,axes=(0,1,2)),
                              k_s,
                              integrator,
                             ) / k_s_squared[...,_np.newaxis]

    return _np.fft.irfftn(displacement,axes=(0,1,2),s=F.shape[:3])


def cell_displacement_avg(size,F):
    """
    Cell center displacement field from average part of the deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.

    """
    F_avg = _np.average(F,axis=(0,1,2))
    return _np.einsum('ml,ijkl->ijkm',F_avg - _np.eye(3),cell_coord0(F.shape[:3],size))


def cell_displacement(size,F):
    """
    Cell center displacement field from deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.

    """
    return cell_displacement_avg(size,F) + cell_displacement_fluct(size,F)


def cell_coord(size,F,origin=_np.zeros(3)):
    """
    Cell center positions.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.
    origin : numpy.ndarray of shape (3), optional
        physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    """
    return cell_coord0(F.shape[:3],size,origin) + cell_displacement(size,F)


def cell_coord0_gridSizeOrigin(coord0,ordered=True):
    """
    Return grid 'DNA', i.e. cells, size, and origin from 1D array of cell positions.

    Parameters
    ----------
    coord0 : numpy.ndarray of shape (:,3)
        undeformed cell coordinates.
    ordered : bool, optional
        expect coord0 data to be ordered (x fast, z slow).

    """
    coords    = [_np.unique(coord0[:,i]) for i in range(3)]
    mincorner = _np.array(list(map(min,coords)))
    maxcorner = _np.array(list(map(max,coords)))
    cells     = _np.array(list(map(len,coords)),'i')
    size      = cells/_np.maximum(cells-1,1) * (maxcorner-mincorner)
    delta     = size/cells
    origin    = mincorner - delta*.5

    # 1D/2D: size/origin combination undefined, set origin to 0.0
    size  [_np.where(cells==1)] = origin[_np.where(cells==1)]*2.
    origin[_np.where(cells==1)] = 0.0

    if cells.prod() != len(coord0):
        raise ValueError('Data count {len(coord0)} does not match cells {cells}.')

    start = origin + delta*.5
    end   = origin - delta*.5 + size

    atol = _np.max(size)*5e-2
    if not (_np.allclose(coords[0],_np.linspace(start[0],end[0],cells[0]),atol=atol) and \
            _np.allclose(coords[1],_np.linspace(start[1],end[1],cells[1]),atol=atol) and \
            _np.allclose(coords[2],_np.linspace(start[2],end[2],cells[2]),atol=atol)):
        raise ValueError('Regular cells spacing violated.')

    if ordered and not _np.allclose(coord0.reshape(tuple(cells)+(3,),order='F'),cell_coord0(cells,size,origin),atol=atol):
        raise ValueError('Input data is not ordered (x fast, z slow).')

    return (cells,size,origin)


def coord0_check(coord0):
    """
    Check whether coordinates lie on a regular grid.

    Parameters
    ----------
    coord0 : numpy.ndarray
        array of undeformed cell coordinates.

    """
    cell_coord0_gridSizeOrigin(coord0,ordered=True)


def node_coord0(cells,size,origin=_np.zeros(3)):
    """
    Nodal positions (undeformed).

    Parameters
    ----------
    cells : numpy.ndarray of shape (3)
        number of cells.
    size : numpy.ndarray of shape (3)
        physical size of the periodic field.
    origin : numpy.ndarray of shape (3), optional
        physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    """
    return _np.stack(_np.meshgrid(_np.linspace(origin[0],size[0]+origin[0],cells[0]+1),
                                  _np.linspace(origin[1],size[1]+origin[1],cells[1]+1),
                                  _np.linspace(origin[2],size[2]+origin[2],cells[2]+1),indexing = 'ij'),
                     axis = -1)


def node_displacement_fluct(size,F):
    """
    Nodal displacement field from fluctuation part of the deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.

    """
    return cell_2_node(cell_displacement_fluct(size,F))


def node_displacement_avg(size,F):
    """
    Nodal displacement field from average part of the deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.

    """
    F_avg = _np.average(F,axis=(0,1,2))
    return _np.einsum('ml,ijkl->ijkm',F_avg - _np.eye(3),node_coord0(F.shape[:3],size))


def node_displacement(size,F):
    """
    Nodal displacement field from deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.

    """
    return node_displacement_avg(size,F) + node_displacement_fluct(size,F)


def node_coord(size,F,origin=_np.zeros(3)):
    """
    Nodal positions.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.
    origin : numpy.ndarray of shape (3), optional
        physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    """
    return node_coord0(F.shape[:3],size,origin) + node_displacement(size,F)


def cell_2_node(cell_data):
    """Interpolate periodic cell data to nodal data."""
    n = (  cell_data + _np.roll(cell_data,1,(0,1,2))
         + _np.roll(cell_data,1,(0,))  + _np.roll(cell_data,1,(1,))  + _np.roll(cell_data,1,(2,))
         + _np.roll(cell_data,1,(0,1)) + _np.roll(cell_data,1,(1,2)) + _np.roll(cell_data,1,(2,0)))*0.125

    return _np.pad(n,((0,1),(0,1),(0,1))+((0,0),)*len(cell_data.shape[3:]),mode='wrap')


def node_2_cell(node_data):
    """Interpolate periodic nodal data to cell data."""
    c = (  node_data + _np.roll(node_data,1,(0,1,2))
         + _np.roll(node_data,1,(0,))  + _np.roll(node_data,1,(1,))  + _np.roll(node_data,1,(2,))
         + _np.roll(node_data,1,(0,1)) + _np.roll(node_data,1,(1,2)) + _np.roll(node_data,1,(2,0)))*0.125

    return c[1:,1:,1:]


def node_coord0_gridSizeOrigin(coord0,ordered=True):
    """
    Return grid 'DNA', i.e. cells, size, and origin from 1D array of nodal positions.

    Parameters
    ----------
    coord0 : numpy.ndarray of shape (:,3)
        undeformed nodal coordinates.
    ordered : bool, optional
        expect coord0 data to be ordered (x fast, z slow).

    """
    coords    = [_np.unique(coord0[:,i]) for i in range(3)]
    mincorner = _np.array(list(map(min,coords)))
    maxcorner = _np.array(list(map(max,coords)))
    cells     = _np.array(list(map(len,coords)),'i') - 1
    size      = maxcorner-mincorner
    origin    = mincorner

    if (cells+1).prod() != len(coord0):
        raise ValueError('Data count {len(coord0)} does not match cells {cells}.')

    atol = _np.max(size)*5e-2
    if not (_np.allclose(coords[0],_np.linspace(mincorner[0],maxcorner[0],cells[0]+1),atol=atol) and \
            _np.allclose(coords[1],_np.linspace(mincorner[1],maxcorner[1],cells[1]+1),atol=atol) and \
            _np.allclose(coords[2],_np.linspace(mincorner[2],maxcorner[2],cells[2]+1),atol=atol)):
        raise ValueError('Regular cells spacing violated.')

    if ordered and not _np.allclose(coord0.reshape(tuple(cells+1)+(3,),order='F'),node_coord0(cells,size,origin),atol=atol):
        raise ValueError('Input data is not ordered (x fast, z slow).')

    return (cells,size,origin)


def regrid(size,F,cells_new):
    """
    Return mapping from coordinates in deformed configuration to a regular cells.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        Physical size.
    F : numpy.ndarray of shape (:,:,:,3,3)
        Deformation gradient field.
    cells_new : numpy.ndarray of shape (3)
        New cells for undeformed coordinates.

    """
    c = cell_coord0(F.shape[:3],size) \
      + cell_displacement_avg(size,F) \
      + cell_displacement_fluct(size,F)

    outer = _np.dot(_np.average(F,axis=(0,1,2)),size)
    for d in range(3):
        c[_np.where(c[:,:,:,d]<0)]        += outer[d]
        c[_np.where(c[:,:,:,d]>outer[d])] -= outer[d]

    tree = _spatial.cKDTree(c.reshape(-1,3),boxsize=outer)
    return tree.query(cell_coord0(cells_new,outer))[1].flatten()
