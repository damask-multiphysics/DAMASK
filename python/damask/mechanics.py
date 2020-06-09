import numpy as _np

def Cauchy(P,F):
    """
    Return Cauchy stress calculated from first Piola-Kirchhoff stress and deformation gradient.

    Resulting tensor is symmetrized as the Cauchy stress needs to be symmetric.

    Parameters
    ----------
    F : numpy.ndarray of shape (...,3,3)
        Deformation gradient.
    P : numpy.ndarray of shape (...,3,3)
        First Piola-Kirchhoff stress.

    """
    sigma = _np.einsum('...,...ij,...kj->...ik',1.0/_np.linalg.det(F),P,F)
    return symmetric(sigma)


def deviatoric_part(T):
    """
    Return deviatoric part of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the deviatoric part is computed.

    """
    return T - _np.einsum('...ij,...->...ij',_np.eye(3),spherical_part(T))


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


def left_stretch(T):
    """
    Return the left stretch of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the left stretch is computed.

    """
    return _polar_decomposition(T,'V')[0]


def maximum_shear(T_sym):
    """
    Return the maximum shear component of a symmetric tensor.

    Parameters
    ----------
    T_sym : numpy.ndarray of shape (...,3,3)
        Symmetric tensor of which the maximum shear is computed.

    """
    w = eigenvalues(T_sym)
    return (w[...,0] - w[...,2])*0.5


def Mises_strain(epsilon):
    """
    Return the Mises equivalent of a strain tensor.

    Parameters
    ----------
    epsilon : numpy.ndarray of shape (...,3,3)
        Symmetric strain tensor of which the von Mises equivalent is computed.

    """
    return _Mises(epsilon,2.0/3.0)


def Mises_stress(sigma):
    """
    Return the Mises equivalent of a stress tensor.

    Parameters
    ----------
    sigma : numpy.ndarray of shape (...,3,3)
        Symmetric stress tensor of which the von Mises equivalent is computed.

    """
    return _Mises(sigma,3.0/2.0)


def PK2(P,F):
    """
    Calculate second Piola-Kirchhoff stress from first Piola-Kirchhoff stress and deformation gradient.

    Parameters
    ----------
    P : numpy.ndarray of shape (...,3,3)
        First Piola-Kirchhoff stress.
    F : numpy.ndarray of shape (...,3,3)
        Deformation gradient.

    """
    S = _np.einsum('...jk,...kl->...jl',_np.linalg.inv(F),P)
    return symmetric(S)


def right_stretch(T):
    """
    Return the right stretch of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the right stretch is computed.

    """
    return _polar_decomposition(T,'U')[0]


def rotational_part(T):
    """
    Return the rotational part of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the rotational part is computed.

    """
    return _polar_decomposition(T,'R')[0]


def spherical_part(T,tensor=False):
    """
    Return spherical (hydrostatic) part of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the hydrostatic part is computed.
    tensor : bool, optional
        Map spherical part onto identity tensor. Default is false

    """
    sph = _np.trace(T,axis2=-2,axis1=-1)/3.0
    return _np.einsum('...jk,...->...jk',_np.eye(3),sph) if tensor else sph


def strain_tensor(F,t,m):
    """
    Return strain tensor calculated from deformation gradient.

    For details refer to https://en.wikipedia.org/wiki/Finite_strain_theory and
    https://de.wikipedia.org/wiki/Verzerrungstensor

    Parameters
    ----------
    F : numpy.ndarray of shape (...,3,3)
        Deformation gradient.
    t : {‘V’, ‘U’}
        Type of the polar decomposition, ‘V’ for left stretch tensor and ‘U’ for right stretch tensor.
    m : float
        Order of the strain.

    """
    if   t == 'V':
        B   = _np.matmul(F,transpose(F))
        w,n = _np.linalg.eigh(B)
    elif t == 'U':
        C   = _np.matmul(transpose(F),F)
        w,n = _np.linalg.eigh(C)

    if   m > 0.0:
        eps = 1.0/(2.0*abs(m)) * (+ _np.matmul(n,_np.einsum('...j,...kj->...jk',w**m,n))
                                  - _np.einsum('...jk->...jk',_np.eye(3)))
        
    elif m < 0.0:
        eps = 1.0/(2.0*abs(m)) * (- _np.matmul(n,_np.einsum('...j,...kj->...jk',w**m,n))
                                  + _np.einsum('...jk->...jk',_np.eye(3)))
    else:
        eps = _np.matmul(n,_np.einsum('...j,...kj->...jk',0.5*_np.log(w),n))

    return eps


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


def _polar_decomposition(T,requested):
    """
    Singular value decomposition.

    Parameters
    ----------
    T : numpy.ndarray of shape (...,3,3)
        Tensor of which the singular values are computed.
    requested : iterable of str
        Requested outputs: ‘R’ for the rotation tensor,
        ‘V’ for left stretch tensor and ‘U’ for right stretch tensor.

    """
    u, s, vh = _np.linalg.svd(T)
    R = _np.einsum('...ij,...jk->...ik',u,vh)

    output = []
    if 'R' in requested:
        output.append(R)
    if 'V' in requested:
        output.append(_np.einsum('...ij,...kj->...ik',T,R))
    if 'U' in requested:
        output.append(_np.einsum('...ji,...jk->...ik',R,T))

    return tuple(output)


def _Mises(T_sym,s):
    """
    Base equation for Mises equivalent of a stres or strain tensor.

    Parameters
    ----------
    T_sym : numpy.ndarray of shape (...,3,3)
        Symmetric tensor of which the von Mises equivalent is computed.
    s : float
        Scaling factor (2/3 for strain, 3/2 for stress).

    """
    d = deviatoric_part(T_sym)
    return _np.sqrt(s*_np.einsum('...jk->...',d**2.0))
