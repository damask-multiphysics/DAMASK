"""
Finite-strain continuum mechanics.

All routines operate on numpy.ndarrays of shape (...,3,3).

"""

from typing import Sequence as _Sequence, Union as _Union #, Literal as _Literal

import numpy as _np

from . import tensor as _tensor
from . import _rotation


def deformation_Cauchy_Green_left(F: _np.ndarray) -> _np.ndarray:
    r"""
    Calculate left Cauchy-Green deformation tensor (Finger deformation tensor).

    Parameters
    ----------
    F : numpy.ndarray, shape (...,3,3)
        Deformation gradient.

    Returns
    -------
    B : numpy.ndarray, shape (...,3,3)
        Left Cauchy-Green deformation tensor.

    Notes
    -----
    .. math::

       \vb{B} = \vb{F} \vb{F}^\text{T}

    """
    return _np.matmul(F,_tensor.transpose(F))


def deformation_Cauchy_Green_right(F: _np.ndarray) -> _np.ndarray:
    r"""
    Calculate right Cauchy-Green deformation tensor.

    Parameters
    ----------
    F : numpy.ndarray, shape (...,3,3)
        Deformation gradient.

    Returns
    -------
    C : numpy.ndarray, shape (...,3,3)
        Right Cauchy-Green deformation tensor.

    Notes
    -----
    .. math::

       \vb{C} = \vb{F}^\text{T} \vb{F}

    """
    return _np.matmul(_tensor.transpose(F),F)


def equivalent_strain_Mises(epsilon: _np.ndarray) -> _np.ndarray:
    r"""
    Calculate the Mises equivalent of a strain tensor.

    Parameters
    ----------
    epsilon : numpy.ndarray, shape (...,3,3)
        Symmetric strain tensor of which the von Mises equivalent is computed.

    Returns
    -------
    epsilon_vM : numpy.ndarray, shape (...)
        Von Mises equivalent strain of epsilon.

    Notes
    -----
    The von Mises equivalent of a strain tensor is defined as:

    .. math::

       \epsilon_\text{vM} = \sqrt{\frac{2}{3}\,\epsilon^\prime_{ij} \epsilon^\prime_{ij}}

    where :math:`\vb*{\epsilon}^\prime` is the deviatoric part
    of the strain tensor.

    """
    return _equivalent_Mises(epsilon,2.0/3.0)


def equivalent_stress_Mises(sigma: _np.ndarray) -> _np.ndarray:
    r"""
    Calculate the Mises equivalent of a stress tensor.

    Parameters
    ----------
    sigma : numpy.ndarray, shape (...,3,3)
        Symmetric stress tensor of which the von Mises equivalent is computed.

    Returns
    -------
    sigma_vM : numpy.ndarray, shape (...)
        Von Mises equivalent stress of sigma.

    Notes
    -----
    The von Mises equivalent of a stress tensor is defined as:

    .. math::

       \sigma_\text{vM} = \sqrt{\frac{3}{2}\,\sigma^\prime_{ij} \sigma^\prime_{ij}}

    where :math:`\vb*{\sigma}^\prime` is the deviatoric part
    of the stress tensor.

    """
    return _equivalent_Mises(sigma,3.0/2.0)


def maximum_shear(T_sym: _np.ndarray) -> _np.ndarray:
    """
    Calculate the maximum shear component of a symmetric tensor.

    Parameters
    ----------
    T_sym : numpy.ndarray, shape (...,3,3)
        Symmetric tensor of which the maximum shear is computed.

    Returns
    -------
    gamma_max : numpy.ndarray, shape (...)
        Maximum shear of T_sym.

    """
    w = _tensor.eigenvalues(T_sym)
    return (w[...,0] - w[...,2])*0.5


def rotation(T: _np.ndarray) -> _rotation.Rotation:
    r"""
    Calculate the rotational part of a tensor.

    Parameters
    ----------
    T : numpy.ndarray, shape (...,3,3)
        Tensor of which the rotational part is computed.

    Returns
    -------
    R : damask.Rotation, shape (...)
        Rotational part of the vector.

    Notes
    -----
    The rotational part is calculated from the polar decomposition:

    .. math::

       \vb{R} = \vb{T} \vb{U}^{-1} = \vb{V}^{-1} \vb{T}

    where :math:`\vb{V}` and :math:`\vb{U}` are the left
    and right stretch tensor, respectively.

    """
    return _rotation.Rotation.from_matrix(_polar_decomposition(T,'R')[0])


def strain(F: _np.ndarray,
           #t: _Literal['V', 'U'],   should work, but rejected by SC
           t: str,
           m: float) -> _np.ndarray:
    r"""
    Calculate strain tensor (Seth–Hill family).

    Parameters
    ----------
    F : numpy.ndarray, shape (...,3,3)
        Deformation gradient.
    t : {'V', 'U'}
        Type of the polar decomposition, 'V' for left stretch tensor
        or 'U' for right stretch tensor.
    m : float
        Order of the strain.

    Returns
    -------
    epsilon : numpy.ndarray, shape (...,3,3)
        Strain of F.

    Notes
    -----
    The strain is defined as:

    .. math::

        \vb*{\epsilon}_V^{(m)} = \frac{1}{2m} (\vb{V}^{2m} - \vb{I}) \\\\
        \vb*{\epsilon}_U^{(m)} = \frac{1}{2m} (\vb{U}^{2m} - \vb{I})

    References
    ----------
    | https://en.wikipedia.org/wiki/Finite_strain_theory
    | https://de.wikipedia.org/wiki/Verzerrungstensor

    """
    if t not in ['V', 'U']: raise ValueError('polar decomposition type not in {V, U}')
    w,n = _np.linalg.eigh(deformation_Cauchy_Green_left(F) if t=='V' else deformation_Cauchy_Green_right(F))
    return    0.5   *  _np.einsum('...j,...kj,...lj',_np.log(w),n,n) if m == 0.0 \
        else  0.5/m * (_np.einsum('...j,...kj,...lj',      w**m,n,n) - _np.eye(3))


def stress_Cauchy(P: _np.ndarray,
                  F: _np.ndarray) -> _np.ndarray:
    """
    Calculate the Cauchy stress (true stress).

    Resulting tensor is symmetrized as the Cauchy stress needs to be symmetric.

    Parameters
    ----------
    P : numpy.ndarray, shape (...,3,3)
        First Piola-Kirchhoff stress.
    F : numpy.ndarray, shape (...,3,3)
        Deformation gradient.

    Returns
    -------
    sigma : numpy.ndarray, shape (...,3,3)
        Cauchy stress.

    """
    return _tensor.symmetric(_np.einsum('...,...ij,...kj',1.0/_np.linalg.det(F),P,F))


def stress_second_Piola_Kirchhoff(P: _np.ndarray,
                                  F: _np.ndarray) -> _np.ndarray:
    """
    Calculate the second Piola-Kirchhoff stress.

    Resulting tensor is symmetrized as the second Piola-Kirchhoff stress
    needs to be symmetric.

    Parameters
    ----------
    P : numpy.ndarray, shape (...,3,3)
        First Piola-Kirchhoff stress.
    F : numpy.ndarray, shape (...,3,3)
        Deformation gradient.

    Returns
    -------
    S : numpy.ndarray, shape (...,3,3)
        Second Piola-Kirchhoff stress.

    """
    return _tensor.symmetric(_np.einsum('...ij,...jk',_np.linalg.inv(F),P))


def stretch_left(T: _np.ndarray) -> _np.ndarray:
    r"""
    Calculate left stretch of a tensor.

    Parameters
    ----------
    T : numpy.ndarray, shape (...,3,3)
        Tensor of which the left stretch is computed.

    Returns
    -------
    V : numpy.ndarray, shape (...,3,3)
        Left stretch tensor from Polar decomposition of T.

    Notes
    -----
    The left stretch tensor is calculated from the
    polar decomposition:

    .. math::

       \vb{V} = \vb{T} \vb{R}^\text{T}

    where :math:`\vb{R}` is a rotation.

    """
    return _polar_decomposition(T,'V')[0]


def stretch_right(T: _np.ndarray) -> _np.ndarray:
    r"""
    Calculate right stretch of a tensor.

    Parameters
    ----------
    T : numpy.ndarray, shape (...,3,3)
        Tensor of which the right stretch is computed.

    Returns
    -------
    U : numpy.ndarray, shape (...,3,3)
        Left stretch tensor from Polar decomposition of T.

    Notes
    -----
    The right stretch tensor is calculated from the
    polar decomposition:

    .. math::

       \vb{U} = \vb{R}^\text{T} \vb{T}

    where :math:`\vb{R}` is a rotation.

    """
    return _polar_decomposition(T,'U')[0]


def _polar_decomposition(T: _np.ndarray,
                         requested: _Union[str, _Sequence[str]]) -> tuple:
    """
    Perform singular value decomposition.

    Parameters
    ----------
    T : numpy.ndarray, shape (...,3,3)
        Tensor of which the singular values are computed.
    requested : sequence of {'R', 'U', 'V'}
        Requested outputs: 'R' for the rotation tensor,
        'V' for left stretch tensor, and 'U' for right stretch tensor.

    Returns
    -------
    VRU : tuple of numpy.ndarray, shape (...,3,3)
       Requested components of the singular value decomposition.

    """
    if isinstance(requested, str): requested = [requested]

    u, _, vh = _np.linalg.svd(T)
    R = _np.einsum('...ij,...jk',u,vh)

    output = []
    if 'R' in requested:
        output+=[R]
    if 'V' in requested:
        output+=[_np.einsum('...ij,...kj',T,R)]
    if 'U' in requested:
        output+=[_np.einsum('...ji,...jk',R,T)]

    if len(output) == 0:
        raise ValueError('output not in {V, R, U}')

    return tuple(output)


def _equivalent_Mises(T_sym: _np.ndarray,
                      s: float) -> _np.ndarray:
    """
    Base equation for Mises equivalent of a stress or strain tensor.

    Parameters
    ----------
    T_sym : numpy.ndarray, shape (...,3,3)
        Symmetric tensor of which the von Mises equivalent is computed.
    s : float
        Scaling factor (2/3 for strain, 3/2 for stress).

    Returns
    -------
    eq : numpy.ndarray, shape (...)
        Scaled second invariant of the deviatoric part of T_sym.

    """
    d = _tensor.deviatoric(T_sym)
    return _np.sqrt(s*_np.sum(d**2.0,axis=(-1,-2)))
