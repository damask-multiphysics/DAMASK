import numpy as np

def Cauchy(P,F):
    """
    Return Cauchy stress calculated from first Piola-Kirchhoff stress and deformation gradient.

    Resulting tensor is symmetrized as the Cauchy stress needs to be symmetric.

    Parameters
    ----------
    F : numpy.ndarray of shape (:,3,3) or (3,3)
        Deformation gradient.
    P : numpy.ndarray of shape (:,3,3) or (3,3)
        First Piola-Kirchhoff stress.

    """
    if np.shape(F) == np.shape(P) == (3,3):
        sigma = 1.0/np.linalg.det(F) * np.dot(P,F.T)
    else:
        sigma = np.einsum('i,ijk,ilk->ijl',1.0/np.linalg.det(F),P,F)
    return symmetric(sigma)


def deviatoric_part(T):
    """
    Return deviatoric part of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (:,3,3) or (3,3)
        Tensor of which the deviatoric part is computed.

    """
    return T - np.eye(3)*spherical_part(T) if np.shape(T) == (3,3) else \
           T - np.einsum('ijk,i->ijk',np.broadcast_to(np.eye(3),[T.shape[0],3,3]),spherical_part(T))


def eigenvalues(T_sym):
    """
    Return the eigenvalues, i.e. principal components, of a symmetric tensor.

    The eigenvalues are sorted in ascending order, each repeated according to
    its multiplicity.

    Parameters
    ----------
    T_sym : numpy.ndarray of shape (:,3,3) or (3,3)
        Symmetric tensor of which the eigenvalues are computed.

    """
    return np.linalg.eigvalsh(symmetric(T_sym))


def eigenvectors(T_sym,RHS=False):
    """
    Return eigenvectors of a symmetric tensor.

    The eigenvalues are sorted in ascending order of their associated eigenvalues.

    Parameters
    ----------
    T_sym : numpy.ndarray of shape (:,3,3) or (3,3)
        Symmetric tensor of which the eigenvectors are computed.
    RHS: bool, optional
        Enforce right-handed coordinate system. Default is False.

    """
    (u,v) = np.linalg.eigh(symmetric(T_sym))

    if RHS:
        if np.shape(T_sym) == (3,3):
            if np.linalg.det(v) < 0.0: v[:,2] *= -1.0
        else:
            v[np.linalg.det(v) < 0.0,:,2] *= -1.0
    return v


def left_stretch(T):
    """
    Return the left stretch of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (:,3,3) or (3,3)
        Tensor of which the left stretch is computed.

    """
    return _polar_decomposition(T,'V')[0]


def maximum_shear(T_sym):
    """
    Return the maximum shear component of a symmetric tensor.

    Parameters
    ----------
    T_sym : numpy.ndarray of shape (:,3,3) or (3,3)
        Symmetric tensor of which the maximum shear is computed.

    """
    w = eigenvalues(T_sym)
    return (w[0] - w[2])*0.5 if np.shape(T_sym) == (3,3) else \
           (w[:,0] - w[:,2])*0.5


def Mises_strain(epsilon):
    """
    Return the Mises equivalent of a strain tensor.

    Parameters
    ----------
    epsilon : numpy.ndarray of shape (:,3,3) or (3,3)
        Symmetric strain tensor of which the von Mises equivalent is computed.

    """
    return _Mises(epsilon,2.0/3.0)


def Mises_stress(sigma):
    """
    Return the Mises equivalent of a stress tensor.

    Parameters
    ----------
    sigma : numpy.ndarray of shape (:,3,3) or (3,3)
        Symmetric stress tensor of which the von Mises equivalent is computed.

    """
    return _Mises(sigma,3.0/2.0)


def PK2(P,F):
    """
    Calculate second Piola-Kirchhoff stress from first Piola-Kirchhoff stress and deformation gradient.

    Parameters
    ----------
    P : numpy.ndarray of shape (:,3,3) or (3,3)
        First Piola-Kirchhoff stress.
    F : numpy.ndarray of shape (:,3,3) or (3,3)
        Deformation gradient.

    """
    if np.shape(F) == np.shape(P) == (3,3):
        S = np.dot(np.linalg.inv(F),P)
    else:
        S = np.einsum('ijk,ikl->ijl',np.linalg.inv(F),P)
    return symmetric(S)


def right_stretch(T):
    """
    Return the right stretch of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (:,3,3) or (3,3)
        Tensor of which the right stretch is computed.

    """
    return _polar_decomposition(T,'U')[0]


def rotational_part(T):
    """
    Return the rotational part of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (:,3,3) or (3,3)
        Tensor of which the rotational part is computed.

    """
    return _polar_decomposition(T,'R')[0]


def spherical_part(T,tensor=False):
    """
    Return spherical (hydrostatic) part of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (:,3,3) or (3,3)
        Tensor of which the hydrostatic part is computed.
    tensor : bool, optional
        Map spherical part onto identity tensor. Default is false

    """
    if T.shape == (3,3):
        sph = np.trace(T)/3.0
        return sph if not tensor else np.eye(3)*sph
    else:
        sph = np.trace(T,axis1=1,axis2=2)/3.0
        if not tensor:
            return sph
        else:
            return np.einsum('ijk,i->ijk',np.broadcast_to(np.eye(3),(T.shape[0],3,3)),sph)


def strain_tensor(F,t,m):
    """
    Return strain tensor calculated from deformation gradient.

    For details refer to https://en.wikipedia.org/wiki/Finite_strain_theory and
    https://de.wikipedia.org/wiki/Verzerrungstensor

    Parameters
    ----------
    F : numpy.ndarray of shape (:,3,3) or (3,3)
        Deformation gradient.
    t : {‘V’, ‘U’}
        Type of the polar decomposition, ‘V’ for left stretch tensor and ‘U’ for right stretch tensor.
    m : float
        Order of the strain.

    """
    F_ = F.reshape(1,3,3) if F.shape == (3,3) else F
    if   t == 'V':
        B   = np.matmul(F_,transpose(F_))
        w,n = np.linalg.eigh(B)
    elif t == 'U':
        C   = np.matmul(transpose(F_),F_)
        w,n = np.linalg.eigh(C)

    if   m > 0.0:
        eps = 1.0/(2.0*abs(m)) * (+ np.matmul(n,np.einsum('ij,ikj->ijk',w**m,n))
                                  - np.broadcast_to(np.eye(3),[F_.shape[0],3,3]))
    elif m < 0.0:
        eps = 1.0/(2.0*abs(m)) * (- np.matmul(n,np.einsum('ij,ikj->ijk',w**m,n))
                                  + np.broadcast_to(np.eye(3),[F_.shape[0],3,3]))
    else:
        eps = np.matmul(n,np.einsum('ij,ikj->ijk',0.5*np.log(w),n))

    return eps.reshape(3,3) if np.shape(F) == (3,3) else \
           eps


def symmetric(T):
    """
    Return the symmetrized tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (:,3,3) or (3,3)
        Tensor of which the symmetrized values are computed.

    """
    return (T+transpose(T))*0.5


def transpose(T):
    """
    Return the transpose of a tensor.

    Parameters
    ----------
    T : numpy.ndarray of shape (:,3,3) or (3,3)
        Tensor of which the transpose is computed.

    """
    return T.T if np.shape(T) == (3,3) else \
           np.transpose(T,(0,2,1))


def _polar_decomposition(T,requested):
    """
    Singular value decomposition.

    Parameters
    ----------
    T : numpy.ndarray of shape (:,3,3) or (3,3)
        Tensor of which the singular values are computed.
    requested : iterable of str
        Requested outputs: ‘R’ for the rotation tensor,
        ‘V’ for left stretch tensor and ‘U’ for right stretch tensor.

    """
    u, s, vh = np.linalg.svd(T)
    R = np.dot(u,vh) if np.shape(T) == (3,3) else \
        np.einsum('ijk,ikl->ijl',u,vh)

    output = []
    if 'R' in requested:
        output.append(R)
    if 'V' in requested:
        output.append(np.dot(T,R.T) if np.shape(T) == (3,3) else np.einsum('ijk,ilk->ijl',T,R))
    if 'U' in requested:
        output.append(np.dot(R.T,T) if np.shape(T) == (3,3) else np.einsum('ikj,ikl->ijl',R,T))

    return tuple(output)


def _Mises(T_sym,s):
    """
    Base equation for Mises equivalent of a stres or strain tensor.

    Parameters
    ----------
    T_sym : numpy.ndarray of shape (:,3,3) or (3,3)
        Symmetric tensor of which the von Mises equivalent is computed.
    s : float
        Scaling factor (2/3 for strain, 3/2 for stress).

    """
    d = deviatoric_part(T_sym)
    return np.sqrt(s*(np.sum(d**2.0))) if np.shape(T_sym) == (3,3) else \
           np.sqrt(s*np.einsum('ijk->i',d**2.0))
