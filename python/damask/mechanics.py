"""Finite-strain continuum mechanics."""

from . import tensor as _tensor
from . import _rotation

import numpy as _np


def deformation_Cauchy_Green_left(F):
    """
    Calculate left Cauchy-Green deformation tensor (Finger deformation tensor).

    Parameters
    ----------
    F : numpy.ndarray of shape (...,3,3)
        Deformation gradient.

    Returns
    -------
    B : numpy.ndarray of shape (...,3,3)
        Left Cauchy-Green deformation tensor.

    """
    return _np.matmul(F,_tensor.transpose(F))


def deformation_Cauchy_Green_right(F):
    """
    Calculate right Cauchy-Green deformation tensor.

    Parameters
    ----------
    F : numpy.ndarray of shape (...,3,3)
        Deformation gradient.

    Returns
    -------
    C : numpy.ndarray of shape (...,3,3)
        Right Cauchy-Green deformation tensor.

    """
    return _np.matmul(_tensor.transpose(F),F)


def equivalent_strain_Mises(epsilon):
    """
    Calculate the Mises equivalent of a strain tensor.

    Parameters
    ----------
    epsilon : numpy.ndarray of shape (...,3,3)
        Symmetric strain tensor of which the von Mises equivalent is computed.

    Returns
    -------
    epsilon_vM : numpy.ndarray of shape (...)
        Von Mises equivalent strain of epsilon.

    """
    return _equivalent_Mises(epsilon,2.0/3.0)


def equivalent_stress_Mises(sigma):
    """
    Calculate the Mises equivalent of a stress tensor.

    Parameters
    ----------
    sigma : numpy.ndarray of shape (...,3,3)
        Symmetric stress tensor of which the von Mises equivalent is computed.

    Returns
    -------
    sigma_vM : numpy.ndarray of shape (...)
        Von Mises equivalent stress of sigma.

    """
    return _equivalent_Mises(sigma,3.0/2.0)


def maximum_shear(T_sym):
    """
    Calculate the maximum shear component of a symmetric tensor.

    Parameters
    ----------
    T_sym : numpy.ndarray of shape (...,3,3)
        Symmetric tensor of which the maximum shear is computed.

    Returns
    -------
    gamma_max : numpy.ndarray of shape (...)
        Maximum shear of T_sym.

    """
    w = _tensor.eigenvalues(T_sym)
    return (w[...,0] - w[...,2])*0.5


def rotation(T):
    """
    Calculate the rotational part of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the rotational part is computed.

    Returns
    -------
    R : damask.Rotation of shape (...)
        Rotational part of the vector.

    """
    return _rotation.Rotation.from_matrix(_polar_decomposition(T,'R')[0])


def strain(F,t,m):
    """
    Calculate strain tensor (Seth–Hill family).

    For details refer to https://en.wikipedia.org/wiki/Finite_strain_theory and
    https://de.wikipedia.org/wiki/Verzerrungstensor

    Parameters
    ----------
    F : numpy.ndarray of shape (...,3,3)
        Deformation gradient.
    t : {‘V’, ‘U’}
        Type of the polar decomposition, ‘V’ for left stretch tensor
        and ‘U’ for right stretch tensor.
    m : float
        Order of the strain.

    Returns
    -------
    epsilon : numpy.ndarray of shape (...,3,3)
        Strain of F.

    """
    if   t == 'V':
        w,n = _np.linalg.eigh(deformation_Cauchy_Green_left(F))
    elif t == 'U':
        w,n = _np.linalg.eigh(deformation_Cauchy_Green_right(F))

    if   m > 0.0:
        eps = 1.0/(2.0*abs(m)) * (+ _np.einsum('...j,...kj,...lj',w**m,n,n) - _np.eye(3))

    elif m < 0.0:
        eps = 1.0/(2.0*abs(m)) * (- _np.einsum('...j,...kj,...lj',w**m,n,n) + _np.eye(3))
    else:
        eps = _np.einsum('...j,...kj,...lj',0.5*_np.log(w),n,n)

    return eps



def stress_Cauchy(P,F):
    """
    Calculate the Cauchy stress (true stress).

    Resulting tensor is symmetrized as the Cauchy stress needs to be symmetric.

    Parameters
    ----------
    P : numpy.ndarray of shape (...,3,3)
        First Piola-Kirchhoff stress.
    F : numpy.ndarray of shape (...,3,3)
        Deformation gradient.

    Returns
    -------
    sigma : numpy.ndarray of shape (...,3,3)
        Cauchy stress.

    """
    return _tensor.symmetric(_np.einsum('...,...ij,...kj',1.0/_np.linalg.det(F),P,F))


def stress_second_Piola_Kirchhoff(P,F):
    """
    Calculate the second Piola-Kirchhoff stress.

    Resulting tensor is symmetrized as the second Piola-Kirchhoff stress
    needs to be symmetric.

    Parameters
    ----------
    P : numpy.ndarray of shape (...,3,3)
        First Piola-Kirchhoff stress.
    F : numpy.ndarray of shape (...,3,3)
        Deformation gradient.

    Returns
    -------
    S : numpy.ndarray of shape (...,3,3)
        Second Piola-Kirchhoff stress.

    """
    return _tensor.symmetric(_np.einsum('...ij,...jk',_np.linalg.inv(F),P))


def stretch_left(T):
    """
    Calculate left stretch of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the left stretch is computed.

    Returns
    -------
    V : numpy.ndarray of shape (...,3,3)
        Left stretch tensor from Polar decomposition of T.

    """
    return _polar_decomposition(T,'V')[0]


def stretch_right(T):
    """
    Calculate right stretch of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the right stretch is computed.

    Returns
    -------
    U : numpy.ndarray of shape (...,3,3)
        Left stretch tensor from Polar decomposition of T.

    """
    return _polar_decomposition(T,'U')[0]


def _polar_decomposition(T,requested):
    """
    Perform singular value decomposition.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the singular values are computed.
    requested : iterable of str
        Requested outputs: ‘R’ for the rotation tensor,
        ‘V’ for left stretch tensor and ‘U’ for right stretch tensor.

    """
    u, _, vh = _np.linalg.svd(T)
    R = _np.einsum('...ij,...jk',u,vh)

    output = []
    if 'R' in requested:
        output.append(R)
    if 'V' in requested:
        output.append(_np.einsum('...ij,...kj',T,R))
    if 'U' in requested:
        output.append(_np.einsum('...ji,...jk',R,T))

    if len(output) == 0:
        raise ValueError('output needs to be out of V, R, U')

    return tuple(output)


def _equivalent_Mises(T_sym,s):
    """
    Base equation for Mises equivalent of a stress or strain tensor.

    Parameters
    ----------
    T_sym : numpy.ndarray of shape (...,3,3)
        Symmetric tensor of which the von Mises equivalent is computed.
    s : float
        Scaling factor (2/3 for strain, 3/2 for stress).

    """
    d = _tensor.deviatoric(T_sym)
    return _np.sqrt(s*_np.sum(d**2.0,axis=(-1,-2)))
