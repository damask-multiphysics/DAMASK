"""
Tensor mathematics.

All routines operate on numpy.ndarrays of shape (...,3,3).
"""

import numpy as _np


def deviatoric(T: _np.ndarray) -> _np.ndarray:
    r"""
    Calculate deviatoric part of a tensor.

    Parameters
    ----------
    T : numpy.ndarray, shape (...,3,3)
        Tensor of which the deviatoric part is computed.

    Returns
    -------
    T' : numpy.ndarray, shape (...,3,3)
        Deviatoric part of T.

    See Also
    --------
    spherical : Calculate spherical part of a tensor.

    Notes
    -----
    The deviatoric part of a tensor is defined as:

    .. math::

        \vb{T}' = \vb{T} - \vb{I}_\text{p},

    where :math:`\vb{I}_\text{p}` is the spherical
    part of the tensor mapped onto identity.
    """
    return T - spherical(T,tensor=True)


def eigenvalues(T_sym: _np.ndarray) -> _np.ndarray:
    r"""
    Calculate eigenvalues of a symmetric tensor.

    Parameters
    ----------
    T_sym : numpy.ndarray, shape (...,3,3)
        Symmetric tensor of which the eigenvalues are computed.

    Returns
    -------
    lambda : numpy.ndarray, shape (...,3)
        Eigenvalues of T_sym sorted in ascending order, each repeated
        according to its multiplicity.

    See Also
    --------
    eigenvectors : Calculate eigenvectors of a symmetric tensor.

    Notes
    -----
    The eigenvalues :math:`\lambda` are defined from:

    .. math::

        \vb{T}_\text{sym} \vb{v}= \lambda \vb{v},

    where :math:`\vb{v}` are the eigenvectors.
    """
    return _np.linalg.eigvalsh(symmetric(T_sym))


def eigenvectors(T_sym: _np.ndarray,
                 RHS: bool = False) -> _np.ndarray:
    r"""
    Calculate eigenvectors of a symmetric tensor.

    Parameters
    ----------
    T_sym : numpy.ndarray, shape (...,3,3)
        Symmetric tensor of which the eigenvectors are computed.
    RHS : bool, optional
        Enforce right-handed coordinate system. Defaults to False.

    Returns
    -------
    v : numpy.ndarray, shape (...,3,3)
        Eigenvectors of T_sym sorted in ascending order of their
        associated eigenvalues.

    See Also
    --------
    eigenvalues : Calculate eigenvalues of a symmetric tensor.

    Notes
    -----
    The eigenvectors :math:`\vb{v}` are defined from:

    .. math::

        \vb{T}_\text{sym} \vb{v}= \lambda \vb{v},

    where :math:`\lambda` are the eigenvalues.
    """
    _,v = _np.linalg.eigh(symmetric(T_sym))

    if RHS: v[_np.linalg.det(v) < 0.0,:,2] *= -1.0
    return v


def spherical(T: _np.ndarray,
              tensor: bool = True) -> _np.ndarray:
    r"""
    Calculate spherical part of a tensor.

    Parameters
    ----------
    T : numpy.ndarray, shape (...,3,3)
        Tensor of which the spherical part is computed.
    tensor : bool, optional
        Map spherical part onto identity tensor. Defaults to True.

    Returns
    -------
    I_p/p : numpy.ndarray, shape (...,3,3) / shape (...,)
        Spherical part of tensor T. I_p is an isotropic tensor.

    See Also
    --------
    deviatoric : Calculate deviatoric part of a tensor.

    Notes
    -----
    The spherical part of a tensor is defined as:

    .. math::

        p &= \text{trace}( \vb{T} )/3 \\
        \vb{I}_p &= \vb{I} \, p
    """
    sph = _np.trace(T,axis2=-2,axis1=-1)/3.0
    return _np.einsum('...jk,...',_np.eye(3),sph) if tensor else sph


def symmetric(T: _np.ndarray) -> _np.ndarray:
    r"""
    Symmetrize tensor.

    Parameters
    ----------
    T : numpy.ndarray, shape (...,3,3)
        Tensor of which the symmetrized values are computed.

    Returns
    -------
    T_sym : numpy.ndarray, shape (...,3,3)
        Symmetrized tensor T.

    Notes
    -----
    The symmetrization of a second-order tensor is defined as:

    .. math::

        \vb{T}_\text{sym} = \left( \vb{T}+ \vb{T}^\text{T} \right)/2
    """
    return (T+transpose(T))*0.5


def transpose(T: _np.ndarray) -> _np.ndarray:
    r"""
    Transpose tensor.

    Parameters
    ----------
    T : numpy.ndarray, shape (...,3,3)
        Tensor of which the transpose is computed.

    Returns
    -------
    T.T : numpy.ndarray, shape (...,3,3)
        Transpose of tensor T.

    Notes
    -----
    The transpose of a second-order tensor is defined as:

    .. math::

        T_{ij}^\text{T} = T_{ji}
    """
    return _np.swapaxes(T,axis2=-2,axis1=-1)
