import numpy as np

def Cauchy(F,P):
  """
  Calculate Cauchy stress from 1st Piola-Kirchhoff stress and deformation gradient.
    
  Resulting tensor is symmetrized as the Cauchy stress needs to be symmetric.
  
  Parameters
  ----------
  F : numpy.array of shape (x,3,3) or (3,3)
     Deformation gradient.
  P : numpy.array of shape (x,3,3) or (3,3)
     1st Piola-Kirchhoff.
  """
  if np.shape(F) == np.shape(P) == (3,3):
    sigma = 1.0/np.linalg.det(F) * np.dot(F,P)
    return (sigma+sigma.T)*0.5
  else:
    sigma =  np.einsum('i,ijk,ilk->ijl',1.0/np.linalg.det(F),P,F)
    return (sigma + np.transpose(sigma,(0,2,1)))*0.5


def deviatoric_part(x):
  """
  Calculate deviatoric part of a tensor.
    
  Parameters
  ----------
  x : numpy.array of shape (x,3,3) or (3,3)
     Tensor.
  """
  if np.shape(x) == (3,3):
    return x - np.eye(3)*np.trace(x)/3.0
  else:
    return x - np.einsum('ijk,i->ijk',np.broadcast_to(np.eye(3),[x.shape[0],3,3]),np.trace(x,axis1=1,axis2=2)/3.0) 


def spherical_part(x):
  """
  Calculate spherical(hydrostatic) part of a tensor.
  
  A single scalar is returned, i.e. the hydrostatic part is not mapped on the 3rd order identity matrix. 
    
  Parameters
  ----------
  x : numpy.array of shape (x,3,3) or (3,3)
     Tensor.
  """
  if np.shape(x) == (3,3):
    return np.trace(x)/3.0
  else:
    return np.trace(x,axis1=1,axis2=2)/3.0
    
    
def Mises_stress(sigma):
  """
  Calculate the Mises equivalent of a stress tensor.
  
  Parameters
  ----------
  sigma : numpy.array of shape (x,3,3) or (3,3)
     Symmetric stress tensor.
  """
  s = deviatoric_part(sigma)
  if np.shape(sigma) == (3,3):
    return np.sqrt(3.0/2.0*np.trace(s))
  else:
    return np.sqrt(np.einsum('ijk->i',s)*3.0/2.0)
    
    
def Mises_strain(epsilon):
  """
  Calculate the Mises equivalent of a strain tensor.
  
  Parameters
  ----------
  sigma : numpy.array of shape (x,3,3) or (3,3)
     Symmetric strain tensor.
  """
  s = deviatoric_part(epsilon)
  if np.shape(epsilon) == (3,3):
    return np.sqrt(2.0/3.0*np.trace(s))
  else:
    return np.sqrt(2.0/3.0*np.einsum('ijk->i',s))
