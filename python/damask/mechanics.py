import numpy as np

def Cauchy(F,P):
    """
    Return Cauchy stress calculated from 1. Piola-Kirchhoff stress and deformation gradient.
      
    Resulting tensor is symmetrized as the Cauchy stress needs to be symmetric.
    
    Parameters
    ----------
    F : numpy.array of shape (:,3,3) or (3,3)
      Deformation gradient.
    P : numpy.array of shape (:,3,3) or (3,3)
      1. Piola-Kirchhoff stress.

    """
    if np.shape(F) == np.shape(P) == (3,3):
        sigma = 1.0/np.linalg.det(F) * np.dot(P,F.T)
    else:
        sigma = np.einsum('i,ijk,ilk->ijl',1.0/np.linalg.det(F),P,F)
    return symmetric(sigma)


def PK2(F,P):
    """
    Return 2. Piola-Kirchhoff stress calculated from 1. Piola-Kirchhoff stress and deformation gradient.
      
    Parameters
    ----------
    F : numpy.array of shape (:,3,3) or (3,3)
      Deformation gradient.
    P : numpy.array of shape (:,3,3) or (3,3)
      1. Piola-Kirchhoff stress.

    """
    if np.shape(F) == np.shape(P) == (3,3):
        S = np.dot(np.linalg.inv(F),P)
    else:
        S = np.einsum('ijk,ikl->ijl',np.linalg.inv(F),P)
    return S


def strain_tensor(F,t,m):
    """
    Return strain tensor calculated from deformation gradient.
      
    For details refer to https://en.wikipedia.org/wiki/Finite_strain_theory and
    https://de.wikipedia.org/wiki/Verzerrungstensor
    
    Parameters
    ----------
    F : numpy.array of shape (:,3,3) or (3,3)
      Deformation gradient.
    t : {‘V’, ‘U’}
      Type of the polar decomposition, ‘V’ for left stretch tensor and ‘U’ for right stretch tensor.
    m : float
      Order of the strain.

    """
    F_ = F.reshape((1,3,3)) if F.shape == (3,3) else F
    if   t == 'U':
        B   = np.matmul(F_,transpose(F_))
        w,n = np.linalg.eigh(B)
    elif t == 'V':
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
      
    return eps.reshape((3,3)) if np.shape(F) == (3,3) else \
           eps


def deviatoric_part(x):
    """
    Return deviatoric part of a tensor.
      
    Parameters
    ----------
    x : numpy.array of shape (:,3,3) or (3,3)
      Tensor of which the deviatoric part is computed.

    """
    return x - np.eye(3)*spherical_part(x) if np.shape(x) == (3,3) else \
           x - np.einsum('ijk,i->ijk',np.broadcast_to(np.eye(3),[x.shape[0],3,3]),spherical_part(x)) 


def spherical_part(x):
    """
    Return spherical (hydrostatic) part of a tensor.
    
    A single scalar is returned, i.e. the hydrostatic part is not mapped on the 3rd order identity
    matrix. 
      
    Parameters
    ----------
    x : numpy.array of shape (:,3,3) or (3,3)
      Tensor of which the hydrostatic part is computed.

    """
    return np.trace(x)/3.0 if np.shape(x) == (3,3) else \
           np.trace(x,axis1=1,axis2=2)/3.0
    
    
def Mises_stress(sigma):
    """
    Return the Mises equivalent of a stress tensor.
    
    Parameters
    ----------
    sigma : numpy.array of shape (:,3,3) or (3,3)
      Symmetric stress tensor of which the von Mises equivalent is computed.

    """
    s = deviatoric_part(sigma)
    return np.sqrt(3.0/2.0*(np.sum(s**2.0))) if np.shape(sigma) == (3,3) else \
           np.sqrt(3.0/2.0*np.einsum('ijk->i',s**2.0))
    
    
def Mises_strain(epsilon):
    """
    Return the Mises equivalent of a strain tensor.
    
    Parameters
    ----------
    epsilon : numpy.array of shape (:,3,3) or (3,3)
      Symmetric strain tensor of which the von Mises equivalent is computed.

    """
    s = deviatoric_part(epsilon)
    return np.sqrt(2.0/3.0*(np.sum(s**2.0))) if np.shape(epsilon) == (3,3) else \
           np.sqrt(2.0/3.0*np.einsum('ijk->i',s**2.0))


def symmetric(x):
    """
    Return the symmetrized tensor.
    
    Parameters
    ----------
    x : numpy.array of shape (:,3,3) or (3,3)
      Tensor of which the symmetrized values are computed.

    """
    return (x+transpose(x))*0.5


def maximum_shear(x):
    """
    Return the maximum shear component of a symmetric tensor.
    
    Parameters
    ----------
    x : numpy.array of shape (:,3,3) or (3,3)
      Symmetric tensor of which the maximum shear is computed.

    """
    w = np.linalg.eigvalsh(symmetric(x))                                                            # eigenvalues in ascending order
    return (w[2] - w[0])*0.5 if np.shape(x) == (3,3) else \
           (w[:,2] - w[:,0])*0.5
    
    
def principal_components(x):
    """
    Return the principal components of a symmetric tensor.
    
    The principal components (eigenvalues) are sorted in descending order, each repeated according to
    its multiplicity.
    
    Parameters
    ----------
    x : numpy.array of shape (:,3,3) or (3,3)
      Symmetric tensor of which the principal compontents are computed.

    """
    w = np.linalg.eigvalsh(symmetric(x))                                                            # eigenvalues in ascending order
    return w[::-1] if np.shape(x) == (3,3) else \
           w[:,::-1]
    
    
def transpose(x):
    """
    Return the transpose of a tensor.
      
    Parameters
    ----------
    x : numpy.array of shape (:,3,3) or (3,3)
      Tensor of which the transpose is computed.

    """
    return x.T if np.shape(x) == (3,3) else \
           np.transpose(x,(0,2,1))


def rotational_part(x):
    """
    Return the rotational part of a tensor.
      
    Parameters
    ----------
    x : numpy.array of shape (:,3,3) or (3,3)
      Tensor of which the rotational part is computed.

    """
    return __polar_decomposition(x,'R')[0]


def left_stretch(x):
    """
    Return the left stretch of a tensor.
      
    Parameters
    ----------
    x : numpy.array of shape (:,3,3) or (3,3)
      Tensor of which the left stretch is computed.

    """
    return __polar_decomposition(x,'V')[0]
  
  
def right_stretch(x):
    """
    Return the right stretch of a tensor.
      
    Parameters
    ----------
    x : numpy.array of shape (:,3,3) or (3,3)
      Tensor of which the right stretch is computed.

    """
    return __polar_decomposition(x,'U')[0]


def __polar_decomposition(x,requested):
    """
    Singular value decomposition.
      
    Parameters
    ----------
    x : numpy.array of shape (:,3,3) or (3,3)
      Tensor of which the singular values are computed.
    requested : iterable of str
      Requested outputs: ‘R’ for the rotation tensor, 
      ‘V’ for left stretch tensor and ‘U’ for right stretch tensor.

    """
    u, s, vh = np.linalg.svd(x)
    R = np.dot(u,vh) if np.shape(x) == (3,3) else \
        np.einsum('ijk,ikl->ijl',u,vh)
    
    output = []
    if 'R' in requested:
        output.append(R)
    if 'V' in requested:
        output.append(np.dot(x,R.T) if np.shape(x) == (3,3) else np.einsum('ijk,ilk->ijl',x,R))
    if 'U' in requested:
        output.append(np.dot(R.T,x) if np.shape(x) == (3,3) else np.einsum('ikj,ikl->ijl',R,x))
    
    return tuple(output)
