"""
Filters for operations on regular grids.

The grids are defined as (x,y,z,...) where x is fastest and z is slowest.
This convention is consistent with the layout in grid vti files.

When converting to/from a plain list (e.g. storage in ASCII table),
the following operations are required for tensorial data:

    - D3 = D1.reshape(cells+(-1,),order='F').reshape(cells+(3,3))
    - D1 = D3.reshape(cells+(-1,)).reshape(-1,9,order='F')

"""

from typing import Tuple as _Tuple

from scipy import spatial as _spatial
import numpy as _np

from ._typehints import FloatSequence as _FloatSequence, IntSequence as _IntSequence


def _ks(size: _FloatSequence,
        cells: _IntSequence,
        first_order: bool = False) -> _np.ndarray:
    """
    Get wave numbers operator.

    Parameters
    ----------
    size : sequence of float, len (3)
        Physical size of the periodic field.
    cells : sequence of int, len (3)
        Number of cells.
    first_order : bool, optional
        Correction for first order derivatives, defaults to False.

    """
    k_sk = _np.where(_np.arange(cells[0])>cells[0]//2,
                     _np.arange(cells[0])-cells[0],_np.arange(cells[0]))/size[0]
    if cells[0]%2 == 0 and first_order: k_sk[cells[0]//2] = 0                                       # Nyquist freq=0 for even cells (Johnson, MIT, 2011)

    k_sj = _np.where(_np.arange(cells[1])>cells[1]//2,
                     _np.arange(cells[1])-cells[1],_np.arange(cells[1]))/size[1]
    if cells[1]%2 == 0 and first_order: k_sj[cells[1]//2] = 0                                       # Nyquist freq=0 for even cells (Johnson, MIT, 2011)

    k_si = _np.arange(cells[2]//2+1)/size[2]

    return _np.stack(_np.meshgrid(k_sk,k_sj,k_si,indexing = 'ij'), axis=-1)


def curl(size: _FloatSequence,
         f: _np.ndarray) -> _np.ndarray:
    u"""
    Calculate curl of a vector or tensor field in Fourier space.

    Parameters
    ----------
    size : sequence of float, len (3)
        Physical size of the periodic field.
    f : numpy.ndarray, shape (:,:,:,3) or (:,:,:,3,3)
        Periodic field of which the curl is calculated.

    Returns
    -------
    ∇ × f : numpy.ndarray, shape (:,:,:,3) or (:,:,:,3,3)
        Curl of f.

    """
    n = _np.prod(f.shape[3:])
    k_s = _ks(size,f.shape[:3],True)

    e = _np.zeros((3, 3, 3))
    e[0, 1, 2] = e[1, 2, 0] = e[2, 0, 1] = +1.0                                                     # Levi-Civita symbol
    e[0, 2, 1] = e[2, 1, 0] = e[1, 0, 2] = -1.0

    f_fourier = _np.fft.rfftn(f,axes=(0,1,2))
    curl_ = (_np.einsum('slm,ijkl,ijkm ->ijks' if n == 3 else
                        'slm,ijkl,ijknm->ijksn',e,k_s,f_fourier)*2.0j*_np.pi)                       # vector 3->3, tensor 3x3->3x3

    return _np.fft.irfftn(curl_,axes=(0,1,2),s=f.shape[:3])


def divergence(size: _FloatSequence,
               f: _np.ndarray) -> _np.ndarray:
    u"""
    Calculate divergence of a vector or tensor field in Fourier space.

    Parameters
    ----------
    size : sequence of float, len (3)
        Physical size of the periodic field.
    f : numpy.ndarray, shape (:,:,:,3) or (:,:,:,3,3)
        Periodic field of which the divergence is calculated.

    Returns
    -------
    ∇ · f : numpy.ndarray, shape (:,:,:,1) or (:,:,:,3)
        Divergence of f.

    """
    n = _np.prod(f.shape[3:])
    k_s = _ks(size,f.shape[:3],True)

    f_fourier = _np.fft.rfftn(f,axes=(0,1,2))
    divergence_ = (_np.einsum('ijkl,ijkl ->ijk' if n == 3 else
                              'ijkm,ijklm->ijkl', k_s,f_fourier)*2.0j*_np.pi)                       # vector 3->1, tensor 3x3->3

    return _np.fft.irfftn(divergence_,axes=(0,1,2),s=f.shape[:3])


def gradient(size: _FloatSequence,
             f: _np.ndarray) -> _np.ndarray:
    u"""
    Calculate gradient of a scalar or vector field in Fourier space.

    Parameters
    ----------
    size : sequence of float, len (3)
        Physical size of the periodic field.
    f : numpy.ndarray, shape (:,:,:,1) or (:,:,:,3)
        Periodic field of which the gradient is calculated.

    Returns
    -------
    ∇ f : numpy.ndarray, shape (:,:,:,3) or (:,:,:,3,3)
        Gradient of f.

    """
    n = _np.prod(f.shape[3:])
    k_s = _ks(size,f.shape[:3],True)

    f_fourier = _np.fft.rfftn(f,axes=(0,1,2))
    gradient_ = (_np.einsum('ijkl,ijkm->ijkm' if n == 1 else
                            'ijkl,ijkm->ijklm',f_fourier,k_s)*2.0j*_np.pi)                          # scalar 1->3, vector 3->3x3

    return _np.fft.irfftn(gradient_,axes=(0,1,2),s=f.shape[:3])


def coordinates0_point(cells: _IntSequence,
                       size: _FloatSequence,
                       origin: _FloatSequence = _np.zeros(3)) -> _np.ndarray:
    """
    Cell center positions (undeformed).

    Parameters
    ----------
    cells : sequence of int, len (3)
        Number of cells.
    size : sequence of float, len (3)
        Physical size of the periodic field.
    origin : sequence of float, len(3), optional
        Physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    Returns
    -------
    x_p_0 : numpy.ndarray, shape (:,:,:,3)
        Undeformed cell center coordinates.

    """
    size_ = _np.array(size,float)
    start = origin         + size_/_np.array(cells,_np.int64)*.5
    end   = origin + size_ - size_/_np.array(cells,_np.int64)*.5

    return _np.stack(_np.meshgrid(_np.linspace(start[0],end[0],cells[0]),
                                  _np.linspace(start[1],end[1],cells[1]),
                                  _np.linspace(start[2],end[2],cells[2]),indexing = 'ij'),
                     axis = -1)


def displacement_fluct_point(size: _FloatSequence,
                             F: _np.ndarray) -> _np.ndarray:
    """
    Cell center displacement field from fluctuation part of the deformation gradient field.

    Parameters
    ----------
    size : sequence of float, len (3)
        Physical size of the periodic field.
    F : numpy.ndarray, shape (:,:,:,3,3)
        Deformation gradient field.

    Returns
    -------
    u_p_fluct : numpy.ndarray, shape (:,:,:,3)
        Fluctuating part of the cell center displacements.

    """
    k_s = _ks(size,F.shape[:3],False)
    k_s_squared = _np.einsum('...l,...l',k_s,k_s)
    k_s_squared[0,0,0] = 1.0

    displacement = -_np.einsum('ijkml,ijkl,l->ijkm',
                              _np.fft.rfftn(F,axes=(0,1,2)),
                              k_s,
                              _np.array([0.5j/_np.pi]*3),
                              ) / k_s_squared[...,_np.newaxis]

    return _np.fft.irfftn(displacement,axes=(0,1,2),s=F.shape[:3])


def displacement_avg_point(size: _FloatSequence,
                           F: _np.ndarray) -> _np.ndarray:
    """
    Cell center displacement field from average part of the deformation gradient field.

    Parameters
    ----------
    size : sequence of float, len (3)
        Physical size of the periodic field.
    F : numpy.ndarray, shape (:,:,:,3,3)
        Deformation gradient field.

    Returns
    -------
    u_p_avg : numpy.ndarray, shape (:,:,:,3)
        Average part of the cell center displacements.

    """
    F_avg = _np.average(F,axis=(0,1,2))
    return _np.einsum('ml,ijkl->ijkm',F_avg - _np.eye(3),coordinates0_point(F.shape[:3],size))


def displacement_point(size: _FloatSequence,
                       F: _np.ndarray) -> _np.ndarray:
    """
    Cell center displacement field from deformation gradient field.

    Parameters
    ----------
    size : sequence of float, len (3)
        Physical size of the periodic field.
    F : numpy.ndarray, shape (:,:,:,3,3)
        Deformation gradient field.

    Returns
    -------
    u_p : numpy.ndarray, shape (:,:,:,3)
        Cell center displacements.

    """
    return displacement_avg_point(size,F) + displacement_fluct_point(size,F)


def coordinates_point(size: _FloatSequence,
                      F: _np.ndarray,
                      origin: _FloatSequence = _np.zeros(3)) -> _np.ndarray:
    """
    Cell center positions.

    Parameters
    ----------
    size : sequence of float, len (3)
        Physical size of the periodic field.
    F : numpy.ndarray, shape (:,:,:,3,3)
        Deformation gradient field.
    origin : sequence of float, len(3), optional
        Physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    Returns
    -------
    x_p : numpy.ndarray, shape (:,:,:,3)
        Cell center coordinates.

    """
    return coordinates0_point(F.shape[:3],size,origin) + displacement_point(size,F)


def cellsSizeOrigin_coordinates0_point(coordinates0: _np.ndarray,
                                       ordered: bool = True) -> _Tuple[_np.ndarray,_np.ndarray,_np.ndarray]:
    """
    Return grid 'DNA', i.e. cells, size, and origin from 1D array of point positions.

    Parameters
    ----------
    coordinates0 : numpy.ndarray, shape (:,3)
        Undeformed cell center coordinates.
    ordered : bool, optional
        Expect coordinates0 data to be ordered (x fast, z slow).
        Defaults to True.

    Returns
    -------
    cells, size, origin : Three numpy.ndarray, each of shape (3)
        Information to reconstruct grid.

    """
    coords    = [_np.unique(coordinates0[:,i]) for i in range(3)]
    mincorner = _np.array(list(map(min,coords)))
    maxcorner = _np.array(list(map(max,coords)))
    cells     = _np.array(list(map(len,coords)),_np.int64)
    size      = cells/_np.maximum(cells-1,1) * (maxcorner-mincorner)
    delta     = size/cells
    origin    = mincorner - delta*.5

    # 1D/2D: size/origin combination undefined, set origin to 0.0
    size  [_np.where(cells == 1)] = origin[_np.where(cells == 1)]*2.
    origin[_np.where(cells == 1)] = 0.0

    if cells.prod() != len(coordinates0):
        raise ValueError(f'data count {len(coordinates0)} does not match cells {cells}')

    start = origin + delta*.5
    end   = origin - delta*.5 + size

    atol = _np.max(size)*5e-2
    if not (_np.allclose(coords[0],_np.linspace(start[0],end[0],cells[0]),atol=atol) and \
            _np.allclose(coords[1],_np.linspace(start[1],end[1],cells[1]),atol=atol) and \
            _np.allclose(coords[2],_np.linspace(start[2],end[2],cells[2]),atol=atol)):
        raise ValueError('non-uniform cell spacing')

    if ordered and not _np.allclose(coordinates0.reshape(tuple(cells)+(3,),order='F'),
                                    coordinates0_point(list(cells),size,origin),atol=atol):
        raise ValueError('input data is not ordered (x fast, z slow)')

    return (cells,size,origin)


def coordinates0_node(cells: _IntSequence,
                      size: _FloatSequence,
                      origin: _FloatSequence = _np.zeros(3)) -> _np.ndarray:
    """
    Nodal positions (undeformed).

    Parameters
    ----------
    cells : sequence of int, len (3)
        Number of cells.
    size : sequence of float, len (3)
        Physical size of the periodic field.
    origin : sequence of float, len(3), optional
        Physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    Returns
    -------
    x_n_0 : numpy.ndarray, shape (:,:,:,3)
        Undeformed nodal coordinates.

    """
    return _np.stack(_np.meshgrid(_np.linspace(origin[0],size[0]+origin[0],cells[0]+1),
                                  _np.linspace(origin[1],size[1]+origin[1],cells[1]+1),
                                  _np.linspace(origin[2],size[2]+origin[2],cells[2]+1),indexing = 'ij'),
                     axis = -1)


def displacement_fluct_node(size: _FloatSequence,
                            F: _np.ndarray) -> _np.ndarray:
    """
    Nodal displacement field from fluctuation part of the deformation gradient field.

    Parameters
    ----------
    size : sequence of float, len (3)
        Physical size of the periodic field.
    F : numpy.ndarray, shape (:,:,:,3,3)
        Deformation gradient field.

    Returns
    -------
    u_n_fluct : numpy.ndarray, shape (:,:,:,3)
        Fluctuating part of the nodal displacements.

    """
    return point_to_node(displacement_fluct_point(size,F))


def displacement_avg_node(size: _FloatSequence,
                          F: _np.ndarray) -> _np.ndarray:
    """
    Nodal displacement field from average part of the deformation gradient field.

    Parameters
    ----------
    size : sequence of float, len (3)
        Physical size of the periodic field.
    F : numpy.ndarray, shape (:,:,:,3,3)
        Deformation gradient field.

    Returns
    -------
    u_n_avg : numpy.ndarray, shape (:,:,:,3)
        Average part of the nodal displacements.

    """
    F_avg = _np.average(F,axis=(0,1,2))
    return _np.einsum('ml,ijkl->ijkm',F_avg - _np.eye(3),coordinates0_node(F.shape[:3],size))


def displacement_node(size: _FloatSequence,
                      F: _np.ndarray) -> _np.ndarray:
    """
    Nodal displacement field from deformation gradient field.

    Parameters
    ----------
    size : sequence of float, len (3)
        Physical size of the periodic field.
    F : numpy.ndarray, shape (:,:,:,3,3)
        Deformation gradient field.

    Returns
    -------
    u_n : numpy.ndarray, shape (:,:,:,3)
        Nodal displacements.

    """
    return displacement_avg_node(size,F) + displacement_fluct_node(size,F)


def coordinates_node(size: _FloatSequence,
                     F: _np.ndarray,
                     origin: _FloatSequence = _np.zeros(3)) -> _np.ndarray:
    """
    Nodal positions.

    Parameters
    ----------
    size : sequence of float, len (3)
        Physical size of the periodic field.
    F : numpy.ndarray, shape (:,:,:,3,3)
        Deformation gradient field.
    origin : sequence of float, len(3), optional
        Physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    Returns
    -------
    x_n : numpy.ndarray, shape (:,:,:,3)
        Nodal coordinates.

    """
    return coordinates0_node(F.shape[:3],size,origin) + displacement_node(size,F)


def cellsSizeOrigin_coordinates0_node(coordinates0: _np.ndarray,
                                      ordered: bool = True) -> _Tuple[_np.ndarray,_np.ndarray,_np.ndarray]:
    """
    Return grid 'DNA', i.e. cells, size, and origin from 1D array of nodal positions.

    Parameters
    ----------
    coordinates0 : numpy.ndarray, shape (:,3)
        Undeformed nodal coordinates.
    ordered : bool, optional
        Expect coordinates0 data to be ordered (x fast, z slow).
        Defaults to True.

    Returns
    -------
    cells, size, origin : Three numpy.ndarray, each of shape (3)
        Information to reconstruct grid.

    """
    coords    = [_np.unique(coordinates0[:,i]) for i in range(3)]
    mincorner = _np.array(list(map(min,coords)))
    maxcorner = _np.array(list(map(max,coords)))
    cells     = _np.array(list(map(len,coords)),_np.int64) - 1
    size      = maxcorner-mincorner
    origin    = mincorner

    if (cells+1).prod() != len(coordinates0):
        raise ValueError(f'data count {len(coordinates0)} does not match cells {cells}')

    atol = _np.max(size)*5e-2
    if not (_np.allclose(coords[0],_np.linspace(mincorner[0],maxcorner[0],cells[0]+1),atol=atol) and \
            _np.allclose(coords[1],_np.linspace(mincorner[1],maxcorner[1],cells[1]+1),atol=atol) and \
            _np.allclose(coords[2],_np.linspace(mincorner[2],maxcorner[2],cells[2]+1),atol=atol)):
        raise ValueError('non-uniform cell spacing')

    if ordered and not _np.allclose(coordinates0.reshape(tuple(cells+1)+(3,),order='F'),
                                    coordinates0_node(list(cells),size,origin),atol=atol):
        raise ValueError('input data is not ordered (x fast, z slow)')

    return (cells,size,origin)


def point_to_node(cell_data: _np.ndarray) -> _np.ndarray:
    """
    Interpolate periodic point data to nodal data.

    Parameters
    ----------
    cell_data : numpy.ndarray, shape (:,:,:,...)
        Data defined on the cell centers of a periodic grid.

    Returns
    -------
    node_data : numpy.ndarray, shape (:,:,:,...)
        Data defined on the nodes of a periodic grid.

    """
    n = (  cell_data + _np.roll(cell_data,1,(0,1,2))
         + _np.roll(cell_data,1,(0,))  + _np.roll(cell_data,1,(1,))  + _np.roll(cell_data,1,(2,))
         + _np.roll(cell_data,1,(0,1)) + _np.roll(cell_data,1,(1,2)) + _np.roll(cell_data,1,(2,0)))*0.125

    return _np.pad(n,((0,1),(0,1),(0,1))+((0,0),)*len(cell_data.shape[3:]),mode='wrap')


def node_to_point(node_data: _np.ndarray) -> _np.ndarray:
    """
    Interpolate periodic nodal data to point data.

    Parameters
    ----------
    node_data : numpy.ndarray, shape (:,:,:,...)
        Data defined on the nodes of a periodic grid.

    Returns
    -------
    cell_data : numpy.ndarray, shape (:,:,:,...)
        Data defined on the cell centers of a periodic grid.

    """
    c = (  node_data + _np.roll(node_data,1,(0,1,2))
         + _np.roll(node_data,1,(0,))  + _np.roll(node_data,1,(1,))  + _np.roll(node_data,1,(2,))
         + _np.roll(node_data,1,(0,1)) + _np.roll(node_data,1,(1,2)) + _np.roll(node_data,1,(2,0)))*0.125

    return c[1:,1:,1:]


def coordinates0_valid(coordinates0: _np.ndarray) -> bool:
    """
    Check whether coordinates form a regular grid.

    Parameters
    ----------
    coordinates0 : numpy.ndarray, shape (:,3)
        Array of undeformed cell coordinates.

    Returns
    -------
    valid : bool
        Whether the coordinates form a regular grid.

    """
    try:
        cellsSizeOrigin_coordinates0_point(coordinates0,ordered=True)
        return True
    except ValueError:
        return False


def unravel_index(idx: _np.ndarray) -> _np.ndarray:
    """
    Convert flat indices to coordinate indices.

    Parameters
    ----------
    idx : numpy.ndarray, shape (:,:,:)
        Grid of flat indices.

    Returns
    -------
    unravelled : numpy.ndarray, shape (:,:,:,3)
        Grid of coordinate indices.

    Examples
    --------
    Unravel a linearly increasing sequence of material indices on a 3 × 2 × 1 grid.

    >>> import numpy as np
    >>> import damask
    >>> seq = np.arange(6).reshape((3,2,1),order='F')
    >>> (coord_idx := damask.grid_filters.unravel_index(seq))
    array([[[[0, 0, 0]],
            [[0, 1, 0]]],
           [[[1, 0, 0]],
            [[1, 1, 0]]],
           [[[2, 0, 0]],
            [[2, 1, 0]]]])
    >>> coord_idx[1,1,0]
    array([1, 1, 0])

    """
    cells = idx.shape
    idx_ = _np.expand_dims(idx,3)
    return _np.block([  idx_ %cells[0],
                       (idx_//cells[0]) %cells[1],
                      ((idx_//cells[0])//cells[1])%cells[2]])

def ravel_index(idx: _np.ndarray) -> _np.ndarray:
    """
    Convert coordinate indices to flat indices.

    Parameters
    ----------
    idx : numpy.ndarray, shape (:,:,:,3)
        Grid of coordinate indices.

    Returns
    -------
    ravelled : numpy.ndarray, shape (:,:,:)
        Grid of flat indices.

    Examples
    --------
    Ravel a reversed sequence of coordinate indices on a 2 × 2 × 1 grid.

    >>> import numpy as np
    >>> import damask
    >>> (rev := np.array([[1,1,0],[0,1,0],[1,0,0],[0,0,0]]).reshape((2,2,1,3)))
    array([[[[1, 1, 0]],
            [[0, 1, 0]]],
           [[[1, 0, 0]],
            [[0, 0, 0]]]])
    >>> (flat_idx := damask.grid_filters.ravel_index(rev))
    array([[[3],
            [2]],
           [[1],
            [0]]])

    """
    cells = idx.shape[:3]
    return   idx[:,:,:,0] \
           + idx[:,:,:,1]*cells[0] \
           + idx[:,:,:,2]*cells[0]*cells[1]


def regrid(size: _FloatSequence,
           F: _np.ndarray,
           cells: _IntSequence) -> _np.ndarray:
    """
    Map a deformed grid A back to a rectilinear grid B.

    The size of grid B is chosen as the average deformed size of grid A.

    Parameters
    ----------
    size : sequence of float, len (3)
        Physical size of grid A.
    F : numpy.ndarray, shape (:,:,:,3,3)
        Deformation gradient field on grid A.
    cells : sequence of int, len (3)
        Cell count along x,y,z of grid B.

    Returns
    -------
    idx : numpy.ndarray of int, shape (cells)
        Flat index of closest point on deformed grid A for each point on grid B.

    """
    box = _np.dot(_np.average(F,axis=(0,1,2)),size)
    c = coordinates_point(size,F)%box
    tree = _spatial.cKDTree(c.reshape((-1,3),order='F'),boxsize=box)
    return tree.query(coordinates0_point(cells,box))[1]
