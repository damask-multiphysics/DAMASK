"""
Filters for operations on regular grids.

Notes
-----
The grids are defined as (x,y,z,...) where x is fastest and z is slowest.
This convention is consistent with the layout in grid vtr files.
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
        Physical size of the periodic field.
    cells : numpy.ndarray of shape (3)
        Number of cells.
    first_order : bool, optional
        Correction for first order derivatives, defaults to False.

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
        Physical size of the periodic field.
    field : numpy.ndarray of shape (:,:,:,3) or (:,:,:,3,3)
        Periodic field of which the curl is calculated.

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
        Physical size of the periodic field.
    field : numpy.ndarray of shape (:,:,:,3) or (:,:,:,3,3)
        Periodic field of which the divergence is calculated.

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
        Physical size of the periodic field.
    field : numpy.ndarray of shape (:,:,:,1) or (:,:,:,3)
        Periodic field of which the gradient is calculated.

    """
    n = _np.prod(field.shape[3:])
    k_s = _ks(size,field.shape[:3],True)

    field_fourier = _np.fft.rfftn(field,axes=(0,1,2))
    grad_ = (_np.einsum('ijkl,ijkm->ijkm', field_fourier,k_s)*2.0j*_np.pi if n == 1 else            # scalar, 1 -> 3
             _np.einsum('ijkl,ijkm->ijklm',field_fourier,k_s)*2.0j*_np.pi)                          # vector, 3 -> 3x3

    return _np.fft.irfftn(grad_,axes=(0,1,2),s=field.shape[:3])


def coordinates0_point(cells,size,origin=_np.zeros(3)):
    """
    Cell center positions (undeformed).

    Parameters
    ----------
    cells : numpy.ndarray of shape (3)
        Number of cells.
    size : numpy.ndarray of shape (3)
        Physical size of the periodic field.
    origin : numpy.ndarray, optional
        Physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    """
    start = origin        + size/cells*.5
    end   = origin + size - size/cells*.5

    return _np.stack(_np.meshgrid(_np.linspace(start[0],end[0],cells[0]),
                                  _np.linspace(start[1],end[1],cells[1]),
                                  _np.linspace(start[2],end[2],cells[2]),indexing = 'ij'),
                     axis = -1)


def displacement_fluct_point(size,F):
    """
    Cell center displacement field from fluctuation part of the deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        Physical size of the periodic field.
    F : numpy.ndarray
        Deformation gradient field.

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


def displacement_avg_point(size,F):
    """
    Cell center displacement field from average part of the deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        Physical size of the periodic field.
    F : numpy.ndarray
        Deformation gradient field.

    """
    F_avg = _np.average(F,axis=(0,1,2))
    return _np.einsum('ml,ijkl->ijkm',F_avg - _np.eye(3),coordinates0_point(F.shape[:3],size))


def displacement_point(size,F):
    """
    Cell center displacement field from deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        Physical size of the periodic field.
    F : numpy.ndarray
        Deformation gradient field.

    """
    return displacement_avg_point(size,F) + displacement_fluct_point(size,F)


def coordinates_point(size,F,origin=_np.zeros(3)):
    """
    Cell center positions.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        Physical size of the periodic field.
    F : numpy.ndarray
        Deformation gradient field.
    origin : numpy.ndarray of shape (3), optional
        Physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    """
    return coordinates0_point(F.shape[:3],size,origin) + displacement_point(size,F)


def cellsSizeOrigin_coordinates0_point(coordinates0,ordered=True):
    """
    Return grid 'DNA', i.e. cells, size, and origin from 1D array of point positions.

    Parameters
    ----------
    coordinates0 : numpy.ndarray of shape (:,3)
        Undeformed cell coordinates.
    ordered : bool, optional
        Expect coordinates0 data to be ordered (x fast, z slow).
        Defaults to True.

    """
    coords    = [_np.unique(coordinates0[:,i]) for i in range(3)]
    mincorner = _np.array(list(map(min,coords)))
    maxcorner = _np.array(list(map(max,coords)))
    cells     = _np.array(list(map(len,coords)),'i')
    size      = cells/_np.maximum(cells-1,1) * (maxcorner-mincorner)
    delta     = size/cells
    origin    = mincorner - delta*.5

    # 1D/2D: size/origin combination undefined, set origin to 0.0
    size  [_np.where(cells==1)] = origin[_np.where(cells==1)]*2.
    origin[_np.where(cells==1)] = 0.0

    if cells.prod() != len(coordinates0):
        raise ValueError(f'Data count {len(coordinates0)} does not match cells {cells}.')

    start = origin + delta*.5
    end   = origin - delta*.5 + size

    atol = _np.max(size)*5e-2
    if not (_np.allclose(coords[0],_np.linspace(start[0],end[0],cells[0]),atol=atol) and \
            _np.allclose(coords[1],_np.linspace(start[1],end[1],cells[1]),atol=atol) and \
            _np.allclose(coords[2],_np.linspace(start[2],end[2],cells[2]),atol=atol)):
        raise ValueError('Regular cell spacing violated.')

    if ordered and not _np.allclose(coordinates0.reshape(tuple(cells)+(3,),order='F'),
                                    coordinates0_point(cells,size,origin),atol=atol):
        raise ValueError('Input data is not ordered (x fast, z slow).')

    return (cells,size,origin)


def coordinates0_check(coordinates0):
    """
    Check whether coordinates lie on a regular grid.

    Parameters
    ----------
    coordinates0 : numpy.ndarray
        Array of undeformed cell coordinates.

    """
    cellsSizeOrigin_coordinates0_point(coordinates0,ordered=True)


def coordinates0_node(cells,size,origin=_np.zeros(3)):
    """
    Nodal positions (undeformed).

    Parameters
    ----------
    cells : numpy.ndarray of shape (3)
        Number of cells.
    size : numpy.ndarray of shape (3)
        Physical size of the periodic field.
    origin : numpy.ndarray of shape (3), optional
        Physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    """
    return _np.stack(_np.meshgrid(_np.linspace(origin[0],size[0]+origin[0],cells[0]+1),
                                  _np.linspace(origin[1],size[1]+origin[1],cells[1]+1),
                                  _np.linspace(origin[2],size[2]+origin[2],cells[2]+1),indexing = 'ij'),
                     axis = -1)


def displacement_fluct_node(size,F):
    """
    Nodal displacement field from fluctuation part of the deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        Physical size of the periodic field.
    F : numpy.ndarray
        Deformation gradient field.

    """
    return point_to_node(displacement_fluct_point(size,F))


def displacement_avg_node(size,F):
    """
    Nodal displacement field from average part of the deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        Physical size of the periodic field.
    F : numpy.ndarray
        Deformation gradient field.

    """
    F_avg = _np.average(F,axis=(0,1,2))
    return _np.einsum('ml,ijkl->ijkm',F_avg - _np.eye(3),coordinates0_node(F.shape[:3],size))


def displacement_node(size,F):
    """
    Nodal displacement field from deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        Physical size of the periodic field.
    F : numpy.ndarray
        Deformation gradient field.

    """
    return displacement_avg_node(size,F) + displacement_fluct_node(size,F)


def coordinates_node(size,F,origin=_np.zeros(3)):
    """
    Nodal positions.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        Physical size of the periodic field.
    F : numpy.ndarray
        Deformation gradient field.
    origin : numpy.ndarray of shape (3), optional
        Physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    """
    return coordinates0_node(F.shape[:3],size,origin) + displacement_node(size,F)


def point_to_node(cell_data):
    """Interpolate periodic point data to nodal data."""
    n = (  cell_data + _np.roll(cell_data,1,(0,1,2))
         + _np.roll(cell_data,1,(0,))  + _np.roll(cell_data,1,(1,))  + _np.roll(cell_data,1,(2,))
         + _np.roll(cell_data,1,(0,1)) + _np.roll(cell_data,1,(1,2)) + _np.roll(cell_data,1,(2,0)))*0.125

    return _np.pad(n,((0,1),(0,1),(0,1))+((0,0),)*len(cell_data.shape[3:]),mode='wrap')


def node_2_point(node_data):
    """Interpolate periodic nodal data to point data."""
    c = (  node_data + _np.roll(node_data,1,(0,1,2))
         + _np.roll(node_data,1,(0,))  + _np.roll(node_data,1,(1,))  + _np.roll(node_data,1,(2,))
         + _np.roll(node_data,1,(0,1)) + _np.roll(node_data,1,(1,2)) + _np.roll(node_data,1,(2,0)))*0.125

    return c[1:,1:,1:]


def cellsSizeOrigin_coordinates0_node(coordinates0,ordered=True):
    """
    Return grid 'DNA', i.e. cells, size, and origin from 1D array of nodal positions.

    Parameters
    ----------
    coordinates0 : numpy.ndarray of shape (:,3)
        Undeformed nodal coordinates.
    ordered : bool, optional
        Expect coordinates0 data to be ordered (x fast, z slow).
        Defaults to True.

    """
    coords    = [_np.unique(coordinates0[:,i]) for i in range(3)]
    mincorner = _np.array(list(map(min,coords)))
    maxcorner = _np.array(list(map(max,coords)))
    cells     = _np.array(list(map(len,coords)),'i') - 1
    size      = maxcorner-mincorner
    origin    = mincorner

    if (cells+1).prod() != len(coordinates0):
        raise ValueError(f'Data count {len(coordinates0)} does not match cells {cells}.')

    atol = _np.max(size)*5e-2
    if not (_np.allclose(coords[0],_np.linspace(mincorner[0],maxcorner[0],cells[0]+1),atol=atol) and \
            _np.allclose(coords[1],_np.linspace(mincorner[1],maxcorner[1],cells[1]+1),atol=atol) and \
            _np.allclose(coords[2],_np.linspace(mincorner[2],maxcorner[2],cells[2]+1),atol=atol)):
        raise ValueError('Regular cell spacing violated.')

    if ordered and not _np.allclose(coordinates0.reshape(tuple(cells+1)+(3,),order='F'),
                                    coordinates0_node(cells,size,origin),atol=atol):
        raise ValueError('Input data is not ordered (x fast, z slow).')

    return (cells,size,origin)


def regrid(size,F,cells):
    """
    Return mapping from coordinates in deformed configuration to a regular grid.

    Parameters
    ----------
    size : numpy.ndarray of shape (3)
        Physical size.
    F : numpy.ndarray of shape (:,:,:,3,3)
        Deformation gradient field.
    cells : numpy.ndarray of shape (3)
        Cell count along x,y,z of remapping grid.

    """
    c = coordinates_point(size,F)

    outer = _np.dot(_np.average(F,axis=(0,1,2)),size)
    for d in range(3):
        c[_np.where(c[:,:,:,d]<0)]        += outer[d]
        c[_np.where(c[:,:,:,d]>outer[d])] -= outer[d]

    tree = _spatial.cKDTree(c.reshape(-1,3),boxsize=outer)
    return tree.query(coordinates0_point(cells,outer))[1].flatten()
