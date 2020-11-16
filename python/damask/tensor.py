"""
Tensor operations.

Notes
-----
This is not a tensor class, but a collection of routines
to operate on numpy.ndarrays of shape (...,3,3).

"""

import numpy as _np


def symmetric(T):
    """
    Symmetrize tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the symmetrized values are computed.

    Returns
    -------
    T_sym : numpy.ndarray of shape (...,3,3)
        Symmetrized tensor T.

    """
    return (T+transpose(T))*0.5


def transpose(T):
    """
    Transpose tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the transpose is computed.

    Returns
    -------
    T.T : numpy.ndarray of shape (...,3,3)
        Transpose of tensor T.

    """
    return _np.swapaxes(T,axis2=-2,axis1=-1)


def eigenvalues(T_sym):
    """
    Eigenvalues, i.e. principal components, of a symmetric tensor.

    Parameters
    ----------
    T_sym : numpy.ndarray of shape (...,3,3)
        Symmetric tensor of which the eigenvalues are computed.

    Returns
    -------
    lambda : numpy.ndarray of shape (...,3)
        Eigenvalues of T_sym sorted in ascending order, each repeated
        according to its multiplicity.

    """
    return _np.linalg.eigvalsh(symmetric(T_sym))


def eigenvectors(T_sym,RHS=False):
    """
    Eigenvectors of a symmetric tensor.

    Parameters
    ----------
    T_sym : numpy.ndarray of shape (...,3,3)
        Symmetric tensor of which the eigenvectors are computed.
    RHS: bool, optional
        Enforce right-handed coordinate system. Defaults to False.

    Returns
    -------
    x : numpy.ndarray of shape (...,3,3)
        Eigenvectors of T_sym sorted in ascending order of their
        associated eigenvalues.

    """
    (u,v) = _np.linalg.eigh(symmetric(T_sym))

    if RHS:
        v[_np.linalg.det(v) < 0.0,:,2] *= -1.0
    return v
