import numpy as np

def Cauchy(F,P):
  if np.shape(F) == np.shape(P) == (3,3):
    sigma = 1.0/np.linalg.det(F) * np.dot(F,P)
    return (sigma+sigma.T)*0.5
  else:
    sigma =  np.einsum('i,ijk,ilk->ijl',1.0/np.linalg.det(F),P,F)
    return (sigma + np.transpose(sigma,(0,2,1)))*0.5


def deviator(x):
  if np.shape(x) == (3,3):
    return x - np.eye(3)*np.trace(x)/3.0
  else:
    return x - np.einsum('ijk,i->ijk',np.broadcast_to(np.eye(3),[x.shape[0],3,3]),np.trace(x,axis1=1,axis2=2)/3.0) 


def spherical(x):
  if np.shape(x) == (3,3):
    return np.trace(x)/3.0
  else:
    return np.trace(x,axis1=1,axis2=2)/3.0
