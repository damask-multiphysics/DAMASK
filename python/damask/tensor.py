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
    Return the symmetrized tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the symmetrized values are computed.

    """
    return (T+transpose(T))*0.5


def transpose(T):
    """
    Return the transpose of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the transpose is computed.

    """
    return _np.swapaxes(T,axis2=-2,axis1=-1)


def eigenvalues(T_sym):
    """
    Return the eigenvalues, i.e. principal components, of a symmetric tensor.

    The eigenvalues are sorted in ascending order, each repeated according to
    its multiplicity.

    Parameters
    ----------
    T_sym : numpy.ndarray of shape (...,3,3)
        Symmetric tensor of which the eigenvalues are computed.

    """
    return _np.linalg.eigvalsh(symmetric(T_sym))


def eigenvectors(T_sym,RHS=False):
    """
    Return eigenvectors of a symmetric tensor.

    The eigenvalues are sorted in ascending order of their associated eigenvalues.

    Parameters
    ----------
    T_sym : numpy.ndarray of shape (...,3,3)
        Symmetric tensor of which the eigenvectors are computed.
    RHS: bool, optional
        Enforce right-handed coordinate system. Default is False.

    """
    (u,v) = _np.linalg.eigh(symmetric(T_sym))

    if RHS:
        v[_np.linalg.det(v) < 0.0,:,2] *= -1.0
    return v
