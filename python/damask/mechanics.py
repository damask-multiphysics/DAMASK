"""Finite-strain continuum mechanics."""

from . import tensor

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
    return _np.matmul(F,tensor.transpose(F))


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
    return _np.matmul(tensor.transpose(F),F)


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
    sigma = _np.einsum('...,...ij,...kj',1.0/_np.linalg.det(F),P,F)
    return tensor.symmetric(sigma)


def deviatoric_part(T):
    """
    Calculate deviatoric part of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the deviatoric part is computed.

    Returns
    -------
    T' : numpy.ndarray of shape (...,3,3)
        Deviatoric part of T.

    """
    return T - _np.einsum('...ij,...',_np.eye(3),spherical_part(T))


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
    w = tensor.eigenvalues(T_sym)
    return (w[...,0] - w[...,2])*0.5


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
    S = _np.einsum('...jk,...kl',_np.linalg.inv(F),P)
    return tensor.symmetric(S)


def rotational_part(T):
    """
    Calculate the rotational part of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the rotational part is computed.

    Returns
    -------
    R : numpy.ndarray of shape (...,3,3)
        Rotational part.

    """
    return _polar_decomposition(T,'R')[0]


def spherical_part(T,tensor=False):
    """
    Calculate spherical (hydrostatic) part of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the hydrostatic part is computed.
    tensor : bool, optional
        Map spherical part onto identity tensor. Defaults to false

    Returns
    -------
    p : numpy.ndarray of shape (...)
        unless tensor == True: shape (...,3,3)
        Spherical part of tensor T, e.g. the hydrostatic part/pressure
        of a stress tensor.

    """
    sph = _np.trace(T,axis2=-2,axis1=-1)/3.0
    return _np.einsum('...jk,...',_np.eye(3),sph) if tensor else sph


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
    d = deviatoric_part(T_sym)
    return _np.sqrt(s*_np.sum(d**2.0,axis=(-1,-2)))


# for compatibility
strain_tensor = strain
