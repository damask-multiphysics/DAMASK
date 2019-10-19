import numpy as np

def Cauchy(F,P):
  """
  Return Cauchy stress calculated from 1. Piola-Kirchhoff stress and deformation gradient.
    
  Resulting tensor is symmetrized as the Cauchy stress needs to be symmetric.
  
  Parameters
  ----------
  F : numpy.array of shape (x,3,3) or (3,3)
    Deformation gradient.
  P : numpy.array of shape (x,3,3) or (3,3)
    1. Piola-Kirchhoff stress.
  """
  if np.shape(F) == np.shape(P) == (3,3):
    sigma = 1.0/np.linalg.det(F) * np.dot(F,P)
  else:
    sigma = np.einsum('i,ijk,ilk->ijl',1.0/np.linalg.det(F),P,F)
  return symmetric(sigma)
  

def strain_tensor(F,t,ord):
  """
  Return strain tensor calculated from deformation gradient.
    
  For details refer to Albrecht Bertram: Elasticity and Plasticity of Large Deformations:
  An Introduction (3rd Edition, 2012), p. 102.
  
  Parameters
  ----------
  F : numpy.array of shape (x,3,3) or (3,3)
    Deformation gradient.
  t : {‘V’, ‘U’}
    Type of the polar decomposition, ‘V’ for right stretch tensor and ‘U’ for left stretch tensor.
  ord : float
    Order of the strain 
  """  
  F_expanded = F if len(F.shape) == 3 else F.reshape(1,3,3)
  
  if   t == 'U':
    B   = np.matmul(F_expanded,transpose(F_expanded))
    U,n = np.linalg.eigh(symmetric(B))
    l   = np.log(U) if ord == 0 else U**ord - np.broadcast_to(np.ones(3),[U.shape[0],3])
  elif t == 'V':
    C   = np.matmul(transpose(F_expanded),F_expanded)
    V,n = np.linalg.eigh(symmetric(C))
    l   = np.log(V) if ord == 0 else np.broadcast_to(np.ones(3),[V.shape[0],3]) - 1.0/V**ord
  
  epsilon = np.matmul(n,np.einsum('ij,ikj->ijk',l,n))
  
  return epsilon.reshape((3,3)) if np.shape(F) == (3,3) else \
         epsilon


def deviatoric_part(x):
  """
  Return deviatoric part of a tensor.
    
  Parameters
  ----------
  x : numpy.array of shape (x,3,3) or (3,3)
    Tensor.
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
  x : numpy.array of shape (x,3,3) or (3,3)
    Tensor.
  """
  
  return np.trace(x)/3.0 if np.shape(x) == (3,3) else \
         np.trace(x,axis1=1,axis2=2)/3.0
    
    
def Mises_stress(sigma):
  """
  Return the Mises equivalent of a stress tensor.
  
  Parameters
  ----------
  sigma : numpy.array of shape (x,3,3) or (3,3)
    Symmetric stress tensor.
  """
  s = deviatoric_part(sigma)
  return np.sqrt(3.0/2.0*np.trace(s)) if np.shape(sigma) == (3,3) else \
         np.sqrt(3.0/2.0*np.einsum('ijk->i',s))
    
    
def Mises_strain(epsilon):
  """
  Return the Mises equivalent of a strain tensor.
  
  Parameters
  ----------
  epsilon : numpy.array of shape (x,3,3) or (3,3)
    Symmetric strain tensor.
  """
  s = deviatoric_part(epsilon)
  return np.sqrt(2.0/3.0*np.trace(s)) if np.shape(epsilon) == (3,3) else \
         np.sqrt(2.0/3.0*np.einsum('ijk->i',s))


def symmetric(x):
  """
  Return the symmetrized tensor.
  
  Parameters
  ----------
  x : numpy.array of shape (x,3,3) or (3,3)
    Tensor.
  """
  return (x+transpose(x))*0.5


def maximum_shear(x):
  """
  Return the maximum shear component of a symmetric tensor.
  
  Parameters
  ----------
  x : numpy.array of shape (x,3,3) or (3,3)
    Symmetric tensor.
  """
  w = np.linalg.eigvalsh(symmetric(x))                                                              # eigenvalues in ascending order
  return (w[2] - w[0])*0.5 if np.shape(epsilon) == (3,3) else \
         (w[:,2] - w[:,0])*0.5
    
    
def principal_components(x):
  """
  Return the principal components of a symmetric tensor.
  
  The principal components (eigenvalues) are sorted in descending order, each repeated according to
  its multiplicity.
  
  Parameters
  ----------
  x : numpy.array of shape (x,3,3) or (3,3)
    Symmetric tensor.
  """
  w = np.linalg.eigvalsh(symmetric(x))                                                              # eigenvalues in ascending order
  return w[::-1] if np.shape(epsilon) == (3,3) else \
         w[:,::-1]
    
    
def transpose(x):
  """
  Return the transpose of a tensor.
    
  Parameters
  ----------
  x : numpy.array of shape (x,3,3) or (3,3)
    Tensor.
  """
  return x.T if np.shape(x) == (3,3) else \
         np.transpose(x,(0,2,1))
