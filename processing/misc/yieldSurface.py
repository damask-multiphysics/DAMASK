#!/usr/bin/python
# -*- coding: UTF-8 no BOM -*-

import threading,time,os,subprocess,shlex,string
import numpy as np
from scipy.linalg import svd
from optparse import OptionParser
import damask
from damask.util import leastsqBound

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = scriptID.split()[1][:-3]

def execute(cmd,streamIn=None,wd='./'):
  '''
    executes a command in given directory and returns stdout and stderr for optional stdin
  '''
  initialPath=os.getcwd()
  os.chdir(wd)
  process = subprocess.Popen(shlex.split(cmd),stdout=subprocess.PIPE,stderr = subprocess.PIPE,stdin=subprocess.PIPE)
  if streamIn != None:
    out,error = process.communicate(streamIn.read())
  else:
    out,error = process.communicate()
  os.chdir(initialPath)

  return out,error

def principalStresses(sigmas):
  '''
    computes principal stresses (i.e. eigenvalues) for a set of Cauchy stresses.
    sorted in descending order.
  '''
  lambdas=np.zeros(0,'d')
  for i in xrange(np.shape(sigmas)[1]):
    eigenvalues = np.linalg.eigvalsh(sym6to33(sigmas[:,i]))
    lambdas = np.append(lambdas,np.sort(eigenvalues)[::-1]) #append eigenvalues in descending order
  lambdas = np.transpose(lambdas.reshape(np.shape(sigmas)[1],3))
  return lambdas

def stressInvariants(lambdas):
  '''
    computes stress invariants (i.e. eigenvalues) for a set of principal Cauchy stresses.
  '''
  Is=np.zeros(0,'d')
  for i in xrange(np.shape(lambdas)[1]):
    I = np.array([lambdas[0,i]+lambdas[1,i]+lambdas[2,i],\
                  lambdas[0,i]*lambdas[1,i]+lambdas[1,i]*lambdas[2,i]+lambdas[2,i]*lambdas[0,i],\
                  lambdas[0,i]*lambdas[1,i]*lambdas[2,i]])
    Is = np.append(Is,I)
  Is = Is.reshape(3,np.shape(lambdas)[1])
  return Is

def formatOutput(n, type='%-14.6f'):
  return ''.join([type for i in xrange(n)])

def sym6to33(sigma6):
  ''' Shape the symmetric stress tensor(6,1) into (3,3) '''
  sigma33 = np.empty((3,3))
  sigma33[0,0] = sigma6[0]; sigma33[1,1] = sigma6[1]; sigma33[2,2] = sigma6[2];
  sigma33[0,1] = sigma6[3]; sigma33[1,0] = sigma6[3]
  sigma33[1,2] = sigma6[4]; sigma33[2,1] = sigma6[4]
  sigma33[2,0] = sigma6[5]; sigma33[0,2] = sigma6[5]
  return sigma33

def array2tuple(array):
  '''transform numpy.array into tuple'''
  try:
    return tuple(array2tuple(i) for i in array)
  except TypeError:
    return array
def get_weight(ndim):
#more to do
  return np.ones(ndim)
# ---------------------------------------------------------------------------------------------
# isotropic yield surfaces
# ---------------------------------------------------------------------------------------------

class Tresca(object):
  '''
    residuum of Tresca yield criterion (eq. 2.26)
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self,sigma0, ydata, sigmas):
    lambdas = principalStresses(sigmas)
    r = np.amax(np.array([abs(lambdas[2,:]-lambdas[1,:]),\
                          abs(lambdas[1,:]-lambdas[0,:]),\
                          abs(lambdas[0,:]-lambdas[2,:])]),0) - sigma0
    return r.ravel()
  def jac(self,sigma0, ydata, sigmas):
    return np.ones(len(ydata)) * (-1.0)

class vonMises(object):
  '''
    residuum of Huber-Mises-Hencky yield criterion (eq. 2.37)
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self, sigma0, ydata, sigmas):
    return HosfordBasis(sigma0, 1.0,1.0,1.0, 2.0, sigmas)
  def jac(self, sigma0, ydata, sigmas):
    return HosfordBasis(sigma0, 1.0,1.0,1.0, 2.0, sigmas, Jac=True, nParas=1)

class Drucker(object):
  '''
    residuum of Drucker yield criterion (eq. 2.41, F = sigma0)
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self, (sigma0, C_D), ydata, sigmas):
    return DruckerBasis(sigma0, C_D, 1.0, sigmas)
  def jac(self, (sigma0, C_D), ydata, sigmas):
    return DruckerBasis(sigma0, C_D, 1.0, sigmas, Jac=True, nParas=2)

class generalDrucker(object):
  '''
    residuum of general Drucker yield criterion (eq. 2.42, F = sigma0)
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self, (sigma0, C_D, p), ydata, sigmas):
    return DruckerBasis(sigma0, C_D, p,  sigmas)
  def jac(self, (sigma0, C_D, p), ydata, sigmas):
    return DruckerBasis(sigma0, C_D, p,  sigmas, Jac=True, nParas=3)

class Hosford(object):
  '''
    residuum of Hershey yield criterion (eq. 2.43, Y = sigma0)
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self, (sigma0, a), ydata, sigmas):
    return HosfordBasis(sigma0, 1.0,1.0,1.0, a, sigmas)
  def jac(self, (sigma0, a), ydata, sigmas):
    return HosfordBasis(sigma0, 1.0,1.0,1.0, a, sigmas, Jac=True, nParas=2)

class Hill1948(object):
  '''
    residuum of Hill 1948 quadratic yield criterion (eq. 2.48)
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self, (F,G,H,L,M,N), ydata, sigmas):
    r = F*(sigmas[1]-sigmas[2])**2.0 + G*(sigmas[2]-sigmas[0])**2.0 + H*(sigmas[0]-sigmas[1])**2.0\
      + 2.0*L*sigmas[4]**2.0 + 2.0*M*sigmas[5]**2.0 + 2.0*N*sigmas[3]**2.0 - 1.0
    return r.ravel()/2.0
  def jac(self, (F,G,H,L,M,N), ydata, sigmas):
    jF=(sigmas[1]-sigmas[2])**2.0; jG=(sigmas[2]-sigmas[0])**2.0; jH=(sigmas[0]-sigmas[1])**2.0
    jL=2.0*sigmas[4]**2.0;         jM=2.0*sigmas[5]**2.0;         jN=2.0*sigmas[3]**2.0
    jaco = []
    for jacv in zip(jF, jG, jH, jL, jM, jN): jaco.append(jacv)
    return np.array(jaco)

class Hill1979(object):
  '''
    residuum of Hill 1979 non-quadratic yield criterion (eq. 2.48)
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self, (f,g,h,a,b,c,m), ydata, sigmas):
    return Hill1979Basis(self.stress0, f,g,h,a,b,c,m, sigmas)
  def jac(self, (f,g,h,a,b,c,m), ydata, sigmas):
    return Hill1979Basis(self.stress0, f,g,h,a,b,c,m, sigmas, Jac=True)

class generalHosford(object):
  '''
    residuum of Hershey yield criterion (eq. 2.104, sigmas = sigma0)
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self, (sigma0, F, G, H, a), ydata, sigmas, nParas=5):
    return HosfordBasis(sigma0, F, G, H, a, sigmas)
  def jac(self, (sigma0, F, G, H, a), ydata, sigmas):
    return HosfordBasis(sigma0, F,G,H, a, sigmas, Jac=True, nParas=5)

class Barlat1991iso(object):
  '''
    residuum of isotropic Barlat 1991 yield criterion (eq. 2.37)
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self, (sigma0, m), ydata, sigmas):
    return Barlat1991Basis(sigma0, 1.0,1.0,1.0,1.0,1.0,1.0, m, sigmas)
  def jac(self, (sigma0, m), ydata, sigmas):
    return Barlat1991Basis(sigma0, 1.0,1.0,1.0,1.0,1.0,1.0, m, sigmas, Jac=True, nParas=2)

class Barlat1991aniso(object):
  '''
    residuum of anisotropic Barlat 1991 yield criterion (eq. 2.37)
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self, (sigma0, a,b,c,f,g,h, m), ydata, sigmas):
    return Barlat1991Basis(sigma0, a,b,c,f,g,h, m, sigmas)
  def jac(self, (sigma0, a,b,c,f,g,h, m), ydata, sigmas):
    return Barlat1991Basis(sigma0, a,b,c,f,g,h, m, sigmas, Jac=True, nParas=8)

class Yld200418p(object):
  '''
    residuum of anisotropic Barlat 1991 yield criterion (eq. 2.37)
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self, (sigma0, c12,c21,c23,c32,c31,c13,c44,c55,c66,
                         d12,d21,d23,d32,d31,d13,d44,d55,d66, m), ydata, sigmas):
    return Yld200418pBasis(sigma0, c12,c21,c23,c32,c31,c13,c44,c55,c66,
                                   d12,d21,d23,d32,d31,d13,d44,d55,d66, m, sigmas)
  def jac(self, (sigma0, c12,c21,c23,c32,c31,c13,c44,c55,c66,
                         d12,d21,d23,d32,d31,d13,d44,d55,d66, m), ydata, sigmas):
    return Yld200418pBasis(sigma0, c12,c21,c23,c32,c31,c13,c44,c55,c66,
                                   d12,d21,d23,d32,d31,d13,d44,d55,d66, m,  sigmas, Jac=True)

class KarafillisBoyce(object):
  '''
    residuum of Karafillis-Boyce yield criterion
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self, (c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,
                 b1, b2, a, alpha), ydata, sigmas):
    return KarafillisBoyceBasis(self.stress0, c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,
                 b1, b2, a, alpha, sigmas)
  def jac(self, (c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,
                 b1, b2, a, alpha), ydata, sigmas):
    return KarafillisBoyceBasis(self.stress0, c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,
                 b1, b2, a, alpha, sigmas, Jac=True)

class BBC2003(object):
  '''
    residuum of anisotropic Barlat 1991 yield criterion (eq. 2.37)
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self, (sigma0, a,b,c, d,e,f,g, k), ydata, sigmas):
    return BBC2003Basis(sigma0, a,b,c, d,e,f,g, k, sigmas)
  def jac(self, (sigma0, a,b,c, d,e,f,g, k), ydata, sigmas):
    return BBC2003Basis(sigma0, a,b,c, d,e,f,g, k, sigmas, Jac=True)

class BBC2005(object):
  '''
    residuum of anisotropic Barlat 1991 yield criterion (eq. 2.37)
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self, (a,b,L, M, N, P, Q, R, k), ydata, sigmas):
    return BBC2005Basis(self.stress0, a,b,L, M, N, P, Q, R, k, sigmas)
  def jac(self, (a,b,L, M, N, P, Q, R, k), ydata, sigmas):
    return BBC2005Basis(self.stress0, a,b,L, M, N, P, Q, R, k, sigmas, Jac=True)

class Cazacu_Barlat2D(object):
  '''
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self, (a1,a2,a3,a4,b1,b2,b3,b4,b5,b10,c), ydata, sigmas):
    return Cazacu_Barlat2DBasis(a1,a2,a3,a4,b1,b2,b3,b4,b5,b10,c, 
                                self.stress0, sigmas)
  def jac(self, (a1,a2,a3,a4,b1,b2,b3,b4,b5,b10,c), ydata, sigmas):
    return Cazacu_Barlat2DBasis(a1,a2,a3,a4,b1,b2,b3,b4,b5,b10,c,  
                                self.stress0, sigmas,Jac=True)

class Cazacu_Barlat3D(object):
  '''
  '''
  def __init__(self, uniaxialStress):
    self.stress0 = uniaxialStress
  def fun(self, (a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,c),ydata, sigmas):
    return Cazacu_Barlat3DBasis(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,c,
                                self.stress0, sigmas)
  def jac(self, (a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,c),ydata, sigmas):
    return Cazacu_Barlat3DBasis(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,c,
                                self.stress0, sigmas,Jac=True)
class Vegter(object):
  '''
    Vegter yield criterion
  '''
  def __init__(self, refPts, refNormals,nspace=11):
    self.refPts, self.refNormals = self._getRefPointsNormals(refPts, refNormals)
    self.hingePts = self._getHingePoints()
    self.nspace = nspace
  def _getRefPointsNormals(self,refPtsQtr,refNormalsQtr):
    if len(refPtsQtr) == 12:
      refPts   = refPtsQtr
      refNormals = refNormalsQtr
    else:
      refPts = np.empty([13,2])
      refNormals = np.empty([13,2])
      refPts[12] = refPtsQtr[0]
      refNormals[12] = refNormalsQtr[0]
      for i in xrange(3):
        refPts[i]   = refPtsQtr[i]
        refPts[i+3] = refPtsQtr[3-i][::-1]
        refPts[i+6] =-refPtsQtr[i]
        refPts[i+9] =-refPtsQtr[3-i][::-1]
        refNormals[i]   = refNormalsQtr[i]
        refNormals[i+3] = refNormalsQtr[3-i][::-1]
        refNormals[i+6] =-refNormalsQtr[i]
        refNormals[i+9] =-refNormalsQtr[3-i][::-1]
    return refPts,refNormals

  def _getHingePoints(self):
    '''
      calculate the hinge point B according to the reference points A,C and the normals n,m
      refPoints  = np.array([[p1_x, p1_y], [p2_x, p2_y]]);
      refNormals = np.array([[n1_x, n1_y], [n2_x, n2_y]])
    '''
    def hingPoint(points, normals):
      A1 = points[0][0];   A2 = points[0][1]
      C1 = points[1][0];   C2 = points[1][1]
      n1 = normals[0][0];  n2 = normals[0][1]
      m1 = normals[1][0];  m2 = normals[1][1]
      B1 = (m2*(n1*A1 + n2*A2) - n2*(m1*C1 + m2*C2))/(n1*m2-m1*n2)
      B2 = (n1*(m1*C1 + m2*C2) - m1*(n1*A1 + n2*A2))/(n1*m2-m1*n2)
      return np.array([B1,B2])
    return np.array([hingPoint(self.refPts[i:i+2],self.refNormals[i:i+2]) for i in xrange(len(self.refPts)-1)])

  def getBezier(self):
    def bezier(R,H):
      b = []
      for mu in np.linspace(0.0,1.0,self.nspace):
        b.append(np.array(R[0]*np.ones_like(mu) + 2.0*mu*(H - R[0]) + mu**2*(R[0]+R[1] - 2.0*H)))
      return b
    return np.array([bezier(self.refPts[i:i+2],self.hingePts[i]) for i in xrange(len(self.refPts)-1)])

def VetgerCriterion(stress,lankford, rhoBi0, theta=0.0):
  '''
    0-pure shear; 1-uniaxial; 2-plane strain; 3-equi-biaxial
  '''
  def getFourierParas(r):
    # get the value after Fourier transformation
    nset = len(r)
    lmatrix = np.empty([nset,nset])
    theta = np.linspace(0.0,np.pi/2,nset)
    for i,th in enumerate(theta):
      lmatrix[i] = np.array([np.cos(2*j*th) for j in xrange(nset)])
    return np.linalg.solve(lmatrix, r)

  nps = len(stress)
  if nps%4 != 0:
    print ('Warning: the number of stress points is uncorrect, stress points of %s are missing in set %i'%(
      ['eq-biaxial, plane strain & uniaxial', 'eq-biaxial & plane strain','eq-biaxial'][nps%4-1],nps/4+1))
  else:
    nset = nps/4
    strsSet = stress.reshape(nset,4,2)
  refPts = np.empty([4,2])

  fouriercoeffs = np.array([np.cos(2.0*i*theta) for i in xrange(nset)])
  for i in xrange(2):
    refPts[3,i] = sum(strsSet[:,3,i])/nset
    for j in xrange(3):
      refPts[j,i] = np.dot(getFourierParas(strsSet[:,j,i]), fouriercoeffs)

  rhoUn = np.dot(getFourierParas(-lankford/(lankford+1)), fouriercoeffs)
  rhoBi = (rhoBi0+1 + (rhoBi0-1)*np.cos(2.0*theta))/(rhoBi0+1 - (rhoBi0-1)*np.cos(2.0*theta))
  nVec = lambda rho : np.array([1.0,rho]/np.sqrt(1.0+rho**2))
  refNormals = np.array([nVec(-1.0),nVec(rhoUn),nVec(0.0),nVec(rhoBi)])

  vegter = Vegter(refPts, refNormals)

def Cazacu_Barlat3DBasis(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,c,
                         sigma0,sigmas, Jac = False):
  '''
    residuum of the 3D Cazacu–Barlat (CZ) yield criterion
  '''
  s11 = sigmas[0]; s22 = sigmas[1]; s33 = sigmas[2]
  s12 = sigmas[3]; s23 = sigmas[4]; s31 = sigmas[5]
  s123, s321 = s11*s22*s33, s12*s23*s31
  s1_2, s2_2, s3_2 = s11**2,   s22**2,   s33**2
  s1_3, s2_3, s3_3 = s11*s1_2, s22*s2_2, s33*s3_2
  s12_2, s23_2, s31_2 = s12**2, s23**2, s31**2
  d12_2, d23_2, d31_2 = (s11-s22)**2, (s22-s33)**2, (s33-s11)**2

  J20 = ( a1*d12_2 + a2*d23_2 + a3*d31_2 )/6.0 + a4*s12_2 + a5*s23_2 + a6*s31_2
  J30 = ( (b1    +b2    )*s1_3  + (b3    +b4    )*s2_3  + ( b1+b4-b2      +  b1+b4-b3     )*s3_3 )/27.0- \
        ( (b1*s22+b2*s33)*s1_2  + (b3*s33+b4*s11)*s2_2  + ((b1+b4-b2)*s11 + (b1+b4-b3)*s22)*s3_2 )/9.0 + \
        ( (b1+b4)*s123/9.0 + b11*s321 )*2.0 - \
        ( ( 2.0*b9 *s22 - b8*s33  - (2.0*b9 -b8)*s11 )*s31_2 + 
          ( 2.0*b10*s33 - b5*s22  - (2.0*b10-b5)*s11 )*s12_2 +
          ( (b6+b7)*s11 - b6*s22  - b7*s33           )*s23_2
        )/3.0
  f0 = (J20**3 - c*J30**2)/18.0
  r  = f0**(1.0/6.0)*(3.0/sigma0)

  if not Jac:
    return (r - 1.0).ravel()
  else:
    drdf = r/f0/108.0
    dj2  =  drdf*3.0*J20**2.0
    dj3  = -drdf*2.0*J30*c
    jc   = -drdf*J30**2

    ja1,ja2,ja3 =  dj2*d12_2/6.0, dj2*d23_2/6.0, dj2*d31_2/6.0
    ja4,ja5,ja6 =  dj2*s12_2,     dj2*s23_2,     dj2*s31_2
    jb1 = dj3*( (s1_3 + 2.0*s3_3)/27.0 - s22*s1_2/9.0 - (s11+s22)*s3_2/9.0 + s123/4.5 )
    jb2 = dj3*( (s1_3 -     s3_3)/27.0 - s33*s1_2/9.0 +  s11     *s3_2/9.0 )
    jb3 = dj3*( (s2_3 -     s3_3)/27.0 - s33*s2_2/9.0 +  s22     *s3_2/9.0 )
    jb4 = dj3*( (s2_3 + 2.0*s3_3)/27.0 - s11*s2_2/9.0 - (s11+s22)*s3_2/9.0 + s123/4.5 )

    jb5, jb10 = dj3*(s22 - s11)*s12_2/3.0,  dj3*(s11 - s33)*s12_2/3.0*2.0
    jb6, jb7  = dj3*(s22 - s11)*s23_2/3.0,  dj3*(s33 - s11)*s23_2/3.0
    jb8, jb9  = dj3*(s33 - s11)*s31_2/3.0,  dj3*(s11 - s22)*s31_2/3.0*2.0
    jb11      = dj3*s321*2.0

    jaco = []
    for jacv in zip(ja1,ja2,ja3,ja4,ja5,ja6,jb1,jb2,jb3,jb4,jb5,jb6,jb7,jb8,jb9,jb10,jb11,jc):
      jaco.append(jacv)
    return np.array(jaco)

def Cazacu_Barlat2DBasis(a1,a2,a3,a4,b1,b2,b3,b4,b5,b10,c, 
                         sigma0,sigmas, Jac = False):
  '''
    residuum of the 2D Cazacu–Barlat (CZ) yield criterion for plain stress
  '''
  s11 = sigmas[0]; s22 = sigmas[1]; s12 = sigmas[3]
  s1_2, s2_2 = s11**2, s22**2
  s1_3, s2_3 = s11*s1_2, s22*s2_2
  s12_2 = s12**2

  J20 = ( a1*(s11-s22)**2 + a2*s2_2 + a3*s1_2 )/6.0 + a4*s12_2
  J30 = ( (b1+b2)*s1_3 + (b3+b4)*s2_3 )/27.0 - ( (b1*s11 + b4*s22)*s11*s22 )/9.0 + \
        (  b5*s22 + (2*b10-b5)*s11 )*s12_2/3.0

  f0 = (J20**3 - c*J30**2)/18.0
  r  = f0**(1.0/6.0)*(3.0/sigma0)

  if not Jac:
    return (r - 1.0).ravel()
  else:
    drdf = r/f0/108.0
    dj2  =  drdf*3.0*J20**2.0
    dj3  = -drdf*2.0*J30*c
    jc   = -drdf*J30**2

    ja1,ja2,ja3,ja4 = dj2*(s11-s22)**2/6.0, dj2*s2_2/6.0, dj2*s1_2/6.0, dj2*s12_2
    jb1, jb2 = s1_3/27.0 - s1_2*s22/9.0,  s1_3/27.0
    jb4, jb3 = s2_3/27.0 - s2_2*s11/9.0,  s2_3/27.0
    jb5, jb10= -s12_2*(s11 - s22)/3.0,    s12_2*s11*2.0/3.0

    jaco = []    
    for jacv in zip(ja1,ja2,ja3,ja4,jb1,jb2,jb3,jb4,jb5,jb10,jc):
      jaco.append(jacv)
    return np.array(jaco)

def DruckerBasis(sigma0, C_D, p, sigmas, Jac=False, nParas=2):
  s11 = sigmas[0]; s22 = sigmas[1]; s33 = sigmas[2]
  s12 = sigmas[3]; s23 = sigmas[4]; s31 = sigmas[5]
  I1  = s11 + s22 + s33
  I2  = s11*s22 + s22*s33 + s33*s11 - s12**2 - s23**2 - s31**2
  I3  = s11*s22*s33 + 2.0*s12*s23*s31 - s12**2*s33 - s23**2*s11 - s31**2*s22
  J2  = I1**2/3.0 - I2
  J3  = I1**3/13.5 - I1*I2/3.0 + I3
  left= J2**(3.0*p) - C_D*J3**(2.0*p); right = 3.0**(0.5)/sigma0
  expo= 1.0/(6.0*p)

  if not Jac:
    return (left**expo*right - 1.0).ravel()
  else:
    jaco = []
    dfdl = expo*left**(expo-1.0)
    js   = -left**expo*right/sigma0
    jC   = -dfdl*J3**(2*p)*right
    if nParas == 2:
      for jacv in zip(js, jC): jaco.append(jacv)
      return np.array(jaco)
    else:
      ln   = lambda x : np.log(x + 1.0e-32)
      dldp = 3.0*J2**(3.0*p)*ln(J2) - 2.0*C_D*J3**(2.0*p)*ln(J3)
      
      jp   = dfdl*dldp*right + (left**expo)*ln(left)*expo/(-p)*right
      for jacv in zip(js, jC, jp): jaco.append(jacv)
      return np.array(jaco)

def Hill1979Basis(sigma0, f,g,h,a,b,c,m, sigmas, Jac=False):

  s1,s2,s3  = principalStresses(sigmas)
  d23 = s2-s3;  d123 = 2.0*s1 - s2 - s3
  d31 = s3-s1;  d231 = 2.0*s2 - s3 - s1
  d12 = s1-s2;  d312 = 2.0*s3 - s1 - s2

  d23s = d23**2; d123s = d123**2
  d31s = d31**2; d231s = d231**2
  d12s = d12**2; d312s = d312**2

  m2 = m/2.0;  mi = 1.0/m
  base = f* d23s**m2 + g* d31s**m2 + h* d12s**m2 + \
         a*d123s**m2 + b*d231s**m2 + c*d312s**m2
  left = base**mi
  r = left/sigma0

  if not Jac:
    return (r-1.0).ravel()
  else:
    ln   = lambda x : np.log(x + 1.0e-32)
    drdb = r/base*mi
    dbdm = ( f* d23s**m2*ln( d23s) + g* d31s**m2*ln( d31s) +  h*d12s**m2*ln( d12s) +
             a*d123s**m2*ln(d123s) + b*d231s**m2*ln(d231s) + c*d312s**m2*ln(d312s) )*0.5
    jf = drdb*d23s**m2;   ja = drdb*d123s**m2
    jg = drdb*d31s**m2;   jb = drdb*d231s**m2
    jh = drdb*d12s**m2;   jc = drdb*d312s**m2
    jm = drdb*dbdm + r*ln(base)*(-mi*mi)

    jaco = []
    for jacv in zip(jf,jg,jh,ja,jb,jc,jm): 
      jaco.append(jacv)
    return np.array(jaco)

def HosfordBasis(sigma0, F,G,H, a, sigmas, Jac=False, nParas=1):
  '''
    residuum of Hershey yield criterion (eq. 2.43, Y = sigma0)
  '''
  lambdas = principalStresses(sigmas)
  diff23  = abs(lambdas[1,:] - lambdas[2,:])
  diff31  = abs(lambdas[2,:] - lambdas[0,:])
  diff12  = abs(lambdas[0,:] - lambdas[1,:])
  base    = F*diff23**a + G*diff31**a + H*diff12**a;   expo = 1.0/a
  left    = base**expo
  right   = 2.0**expo*sigma0

  if not Jac:
    if nParas == 1: return (left - right).ravel()
    else:           return (left/right - 1.0).ravel()
  else:
    ones = np.ones(np.shape(sigmas)[1])
    if nParas > 1:
      ln   = lambda x : np.log(x + 1.0e-32)
      dbda = F*ln(diff23)*diff23**a + G*ln(diff31)*diff31**a + H*ln(diff12)*diff12**a
      deda = -expo*expo
      drda = sigma0*(2.0**expo)*ln(2.0)*deda
      dldb = expo*left/base
      jaco = []

    if   nParas == 1:     # von Mises
      return ones*(-2.0**0.5)
    elif nParas == 2:     # isotropic Hosford
      js = ones*(-2.0**expo)                                              # d[]/dsigma0
      ja = dldb*dbda + left*ln(base)*deda - drda                          # d[]/da
      for jacv in zip(js, ja): 
        jaco.append(jacv)
      return np.array(jaco)
    elif nParas == 5:           # anisotropic Hosford
      js = -left/right/sigma0   #ones*(-2.0**expo)                                               # d[]/dsigma0
      jF = dldb*diff23**a/right
      jG = dldb*diff31**a/right
      jH = dldb*diff12**a/right
      ja =(dldb*dbda + left*ln(base)*deda)/right + left*(-right**(-2))*drda    # d[]/da
      for jacv in zip(js, jF,jG,jH,ja): 
        jaco.append(jacv)
      return np.array(jaco)

def Barlat1991Basis(sigma0, a, b, c, f, g, h, m, sigmas, Jac=False, nParas=2):
  '''
    residuum of Barlat 1997 yield criterion
  '''
  cos  = np.cos; sin = np.sin;  pi = np.pi;    abs = np.abs
  dAda = sigmas[1] - sigmas[2];   A = a*dAda
  dBdb = sigmas[2] - sigmas[0];   B = b*dBdb
  dCdc = sigmas[0] - sigmas[1];   C = c*dCdc
  dFdf = sigmas[4];               F = f*dFdf
  dGdg = sigmas[5];               G = g*dGdg
  dHdh = sigmas[3];               H = h*dHdh

  I2 =   (F*F + G*G  +  H*H)/3.0  + ((A-C)**2+(C-B)**2+(B-A)**2)/54.0
  I3 =   (C-B)*(A-C)*(B-A)/54.0 + F*G*H - \
       ( (C-B)*F*F + (A-C)*G*G  + (B-A)*H*H )/6.0
  theta = np.arccos(I3/I2**1.5)
  phi1 = (2.0*theta + pi)/6.0
  phi2 = (2.0*theta + pi*3.0)/6.0
  phi3 = (2.0*theta + pi*5.0)/6.0
  cos1 = 2.0*cos(phi1);  absc1 = abs(cos1)
  cos2 = 2.0*cos(phi2);  absc2 = abs(cos2)
  cos3 = 2.0*cos(phi3);  absc3 = abs(cos3)
  ratio= np.sqrt(3.0*I2)/sigma0;    expo = 1.0/m
  left = ( absc1**m + absc2**m + absc3**m )/2.0
  leftNorm = left**expo
  r    = ratio*leftNorm - 1.0

  if not Jac:
    return r.ravel()
  else:
    ln   = lambda x : np.log(x + 1.0e-32)
    jaco = []
    dfdl = expo*leftNorm/left
    js = -(r + 1.0)/sigma0
    jm = (r+1.0)*ln(left)*(-expo*expo) + ratio*dfdl*0.5*(
         absc1**m*ln(absc1) + absc2**m*ln(absc2) + absc3**m*ln(absc3) )
    if nParas == 2:
      for jacv in zip(js, jm): jaco.append(jacv)
      return np.array(jaco)
    else:
      dI2da = (2.0*A-B-C)*dAda/27.0
      dI2db = (2.0*B-C-A)*dBdb/27.0
      dI2dc = (2.0*C-A-B)*dCdc/27.0
      dI2df = 2.0*F*dFdf/3.0
      dI2dg = 2.0*G*dGdg/3.0
      dI2dh = 2.0*H*dHdh/3.0
      dI3da = dI2da*(B-C)/2.0 + (H**2 - G**2)*dAda/6.0
      dI3db = dI2db*(C-A)/2.0 + (F**2 - H**2)*dBdb/6.0
      dI3dc = dI2dc*(A-B)/2.0 + (G**2 - F**2)*dCdc/6.0
      dI3df = ( (H*G + (B-C)) * F/3.0 )*dFdf
      dI3dg = ( (F*H + (C-A)) * G/3.0 )*dGdg
      dI3dh = ( (G*F + (A-B)) * H/3.0 )*dHdh

      darccos = -(1.0 - I3**2/I2**3)**(-0.5)
      dthedI2 = darccos*I3*(-1.5)*I2**(-2.5)
      dthedI3 = darccos*I2**(-1.5)
      dc1dthe = -sin(phi1)/3.0
      dc2dthe = -sin(phi2)/3.0
      dc3dthe = -sin(phi3)/3.0
      dfdc  = ratio * dfdl * 0.5 * m
      dfdc1 = dfdc* absc1**(expo-1.0)*np.sign(cos1)
      dfdc2 = dfdc* absc2**(expo-1.0)*np.sign(cos2)
      dfdc3 = dfdc* absc3**(expo-1.0)*np.sign(cos3)
      dfdthe= (dfdc1*dc1dthe + dfdc2*dc2dthe + dfdc2*dc2dthe)*2.0
      dfdI2 = dfdthe*dthedI2;   dfdI3 = dfdthe*dthedI3
      ja = dfdI2*dI2da + dfdI3*dI3da
      jb = dfdI2*dI2db + dfdI3*dI3db
      jc = dfdI2*dI2dc + dfdI3*dI3dc
      jf = dfdI2*dI2df + dfdI3*dI3df
      jg = dfdI2*dI2dg + dfdI3*dI3dg
      jh = dfdI2*dI2dh + dfdI3*dI3dh

      for jacv in zip(js,ja,jb,jc,jf,jg,jh,jm):
        jaco.append(jacv)
      return np.array(jaco)

def BBC2003Basis(sigma0, a,b,c, d,e,f,g, k, sigmas, Jac=False):
  '''
    residuum of the BBC2003 yield criterion for plain stress
  '''
  s11 = sigmas[0]; s22 = sigmas[1]; s12 = sigmas[3]
  k2  = 2.0*k
  M = d+e;  N = e+f;  P = (d-e)/2.0;  Q = (e-f)/2.0; R = g**2
  Gamma =    M*s11 + N*s22
  Psi   = ( (P*s11 + Q*s22)**2 + s12**2*R )**0.5

  l1  = b*Gamma + c*Psi; l2  = b*Gamma - c*Psi; l3  = 2.0*c*Psi
  l1s = l1**2;           l2s = l2**2;           l3s = l3**2
  left = a*l1s**k + a*l2s**k + (1-a)*l3s**k
  sBar = left**(1.0/k2);  r = sBar/sigma0 - 1.0
  if not Jac:
    return r.ravel()
  else:
    temp = (P*s11 + Q*s22)/Psi
    dPsidP = temp*s11;  dPsidQ = temp*s22;  dPsidR = 0.5*s12**2/Psi
    ln   = lambda x : np.log(x + 1.0e-32)
    jaco = []
    expo = 0.5/k;  k1 = k-1.0

    dsBardl = expo*sBar/left/sigma0
    dsBarde = sBar*ln(left);   dedk = expo/(-k)
    dldl1 =    a *k*(l1s**k1)*(2.0*l1)
    dldl2 =    a *k*(l2s**k1)*(2.0*l2)
    dldl3 = (1-a)*k*(l3s**k1)*(2.0*l3)

    dldGama = (dldl1 + dldl2)*b
    dldPsi  = (dldl1 - dldl2 + 2.0*dldl3)*c

    dlda = l1s**k + l2s**k - l3s**k
    dldb = dldl1*Gamma + dldl2*Gamma
    dldc = dldl1*Psi   - dldl2*Psi + dldl3*2.0*Psi
    dldk = a*ln(l1s)*l1s**k + a*ln(l2s)*l2s**k + (1-a)*ln(l3s)*l3s**k

    js = -(r + 1.0)/sigma0
    ja = dsBardl * dlda
    jb = dsBardl * dldb
    jc = dsBardl * dldc
    jd = dsBardl *(dldGama*s11 + dldPsi*dPsidP*0.5)
    je = dsBardl *(dldGama*(s11+s22) + dldPsi*(dPsidP*(-0.5) + dPsidQ*0.5) )
    jf = dsBardl *(dldGama*s22 + dldPsi*dPsidQ*(-0.5))
    jg = dsBardl * dldPsi * dPsidR * 2.0*g
    jk = dsBardl * dldk + dsBarde * dedk

    for jacv in zip(js,ja,jb,jc,jd, je, jf,jg,jk):
      jaco.append(jacv)
    return np.array(jaco)

def BBC2005Basis(sigma0, a,b,L, M, N, P, Q, R, k, sigmas, Jac=False):
  '''
    residuum of the BBC2005 yield criterion for plain stress
  '''
  s11 = sigmas[0]; s22 = sigmas[1]; s12 = sigmas[3]
  k2  = 2.0*k
  Gamma  =    L*s11 + M*s22
  Lambda = ( (N*s11 - P*s22)**2 + s12**2 )**0.5
  Psi    = ( (Q*s11 - R*s22)**2 + s12**2 )**0.5

  l1  = Lambda + Gamma; l2  = Lambda - Gamma; l3  = Lambda + Psi; l4  = Lambda - Psi
  l1s = l1**2;          l2s = l2**2;          l3s = l3**2;        l4s = l4**2
  left = a*l1s**k + a*l2s**k + b*l3s**k + b*l4s**k
  sBar = left**(1.0/k2);  r = sBar/sigma0 - 1.0
  if not Jac:
    return r.ravel()
  else:
    ln   = lambda x : np.log(x + 1.0e-32)
    jaco = []
    expo = 0.5/k;  k1 = k-1.0

    dsBardl = expo*sBar/left/sigma0
    dsBarde = sBar*ln(left);   dedk = expo/(-k)
    dldl1 = a*k*(l1s**k1)*(2.0*l1)
    dldl2 = a*k*(l2s**k1)*(2.0*l2)
    dldl3 = b*k*(l3s**k1)*(2.0*l3)
    dldl4 = b*k*(l4s**k1)*(2.0*l4)
    
    dldLambda = dldl1 + dldl2 + dldl3 + dldl4
    dldGama   = dldl1 - dldl2
    dldPsi    = dldl3 - dldl4
    temp = (N*s11 - P*s22)/Lambda
    dLambdadN = s11*temp; dLambdadP = -s22*temp
    temp = (Q*s11 - R*s22)/Psi
    dPsidQ = s11*temp; dPsidR = -s22*temp
    dldk = a*ln(l1s)*l1s**k + a*ln(l2s)*l2s**k + b*ln(l3s)*l3s**k + b*ln(l4s)*l4s**k

    ja = dsBardl * (l1s**k + l2s**k)
    jb = dsBardl * (l3s**k + l4s**k)
    jL = dsBardl * dldGama*s11
    jM = dsBardl * dldGama*s22
    jN = dsBardl * dldLambda*dLambdadN
    jP = dsBardl * dldLambda*dLambdadP
    jQ = dsBardl * dldPsi*dPsidQ
    jR = dsBardl * dldPsi*dPsidR
    jk = dsBardl * dldk + dsBarde * dedk

    for jacv in zip(ja,jb,jL,jM, jN, jP,jQ,jR,jk):
      jaco.append(jacv)
    return np.array(jaco)

def principalStress(p):
  sin = np.sin;    cos = np.cos
  s11 = p[0];   s22 = p[1];   s33 = p[2]
  s12 = p[3];   s23 = p[4];   s31 = p[5]
  I1  = s11 + s22 + s33
  I2  = s11*s22 + s22*s33 + s33*s11 - s12**2 - s23**2 - s31**2
  I3  = s11*s22*s33 + 2.0*s12*s23*s31 - s12**2*s33 - s23**2*s11 - s31**2*s22

  third = 1.0/3.0
  I1s3I2= (I1**2 - 3.0*I2)**0.5
  numer = 2.0*I1**3 - 9.0*I1*I2 + 27.0*I3
  denom = I1s3I2**(-3.0)
  cs    = 0.5*numer*denom
  phi   = np.arccos(cs)/3.0
  t1    = I1/3.0; t2 = 2.0/3.0*I1s3I2
  S     = np.array( [t1 + t2*cos(phi), t1+t2*cos(phi+np.pi*2.0/3.0), t1+t2*cos(phi+np.pi*4.0/3.0)])
  return S, np.array([I1,I2,I3])

def principalStrs_Der(p, Invariant, s1, s2, s3, s4, s5, s6, Karafillis=False):
  sin = np.sin;    cos = np.cos
  I1 = Invariant[0,:];  I2 = Invariant[1,:];  I3 = Invariant[2,:]

  third = 1.0/3.0
  I1s3I2= (I1**2 - 3.0*I2)**0.5
  numer = 2.0*I1**3 - 9.0*I1*I2 + 27.0*I3
  denom = I1s3I2**(-3.0)
  cs    = 0.5*numer*denom
  phi   = np.arccos(cs)*third

  dphidcs   = -third/np.sqrt(1.0 - cs**2)
  dcsddenom = 0.5*numer*(-1.5)*I1s3I2**(-5.0)
  dcsdI1    = 0.5*(6.0*I1**2 - 9.0*I2)*denom + dcsddenom*(2.0*I1)
  dcsdI2    = 0.5*(          - 9.0*I1)*denom + dcsddenom*(-3.0)
  dcsdI3    = 13.5*denom
  dphidI1, dphidI2, dphidI3 = dphidcs*dcsdI1, dphidcs*dcsdI2, dphidcs*dcsdI3

  dI1s3I2dI1= I1/I1s3I2;  dI1s3I2dI2 = -1.5/I1s3I2
  third2    = 2.0*third;     tcoeff  = third2*I1s3I2
  
  dSidIj = lambda theta : ( tcoeff*(-sin(theta))*dphidI1 + third2*dI1s3I2dI1*cos(theta) + third, 
                            tcoeff*(-sin(theta))*dphidI2 + third2*dI1s3I2dI2*cos(theta), 
                            tcoeff*(-sin(theta))*dphidI3)
  dSdI = np.array([dSidIj(phi),dSidIj(phi+np.pi*2.0/3.0),dSidIj(phi+np.pi*4.0/3.0)]) # i=1,2,3; j=1,2,3
  
# calculate the derivation of principal stress with regards to the anisotropic coefficients  
  one    = np.ones_like(p[0]); zero = np.zeros_like(p[0])
  dIdp = np.array([ [one,    one,    one,    zero,   zero,   zero], 
                    [p[1]+p[2], p[2]+p[0], p[0]+p[1], -2.0*p[3], -2.0*p[4], -2.0*p[5]],
                    [p[1]*p[2]-p[4]**2, p[2]*p[0]-p[5]**2, p[0]*p[1]-p[3]**2, 
                     -2.0*p[3]*p[2]+2.0*p[4]*p[5], -2.0*p[4]*p[0]+2.0*p[5]*p[3], -2.0*p[5]*p[1]+2.0*p[3]*p[4]]
                  ]) 

  if Karafillis:
    dSdp = np.empty_like(dIdp)
    dSdc = np.empty_like(dIdp)
    zero   = np.zeros_like(s1)
    dpdc = np.array([[zero,s2-s3,s3-s2], [s1-s3,zero,s3-s1], [s1-s2,s2-s1,zero]])
    for i in xrange(3):
      for j in xrange(6):
        dSdp[i,j] = dSdI[i,0]*dIdp[0,j]+dSdI[i,1]*dIdp[1,j]+dSdI[i,2]*dIdp[2,j]
      for j in xrange(3):
        dSdc[i,j] = (dSdp[i,0]*dpdc[j,0]+dSdp[i,1]*dpdc[j,1]+dSdp[i,2]*dpdc[j,2])/3.0
      dSdc[i,3:6] = dSdp[i,3]*s4,dSdp[i,4]*s5,dSdp[i,5]*s6
    return dSdc
  else:
    dIdc = np.empty([3,9,len(s1)])
    dSdc = np.empty_like(dIdc)

    for i in xrange(3):
      dIdc[i]=np.array([dIdp[i,0]*(-s2), dIdp[i,1]*(-s1), dIdp[i,1]*(-s3),
                        dIdp[i,2]*(-s2), dIdp[i,2]*(-s1), dIdp[i,0]*(-s3),
                        dIdp[i,3]*  s4,  dIdp[i,4]*  s5,  dIdp[i,5]*  s6 ])
    for i in xrange(3):
      for j in xrange(9):
        dSdc[i,j] = dSdI[i,0]*dIdc[0,j]+dSdI[i,1]*dIdc[1,j]+dSdI[i,2]*dIdc[2,j]
    return dSdc

def Yld200418pBasis(sigma0, c12,c21,c23,c32,c31,c13,c44,c55,c66,
                            d12,d21,d23,d32,d31,d13,d44,d55,d66, m, sigmas, Jac=False):

  sv = (sigmas[0] + sigmas[1] + sigmas[2])/3.0
  s1 = sigmas[0]-sv; s2 = sigmas[1]-sv; s3 = sigmas[2]-sv
  s4 = sigmas[3];    s5 = sigmas[4];    s6 = sigmas[5]
  
  ys = lambda s1,s2,s3,s4,s5,s6,c12,c21,c23,c32,c13,c31,c44,c55,c66: np.array( [
         -c12*s2-c13*s3, -c21*s1-c23*s3, -c31*s1-c32*s2, c44*s4, c55*s5, c66*s6 ])
  p = ys(s1,s2,s3,s4,s5,s6,c12,c21,c23,c32,c13,c31,c44,c55,c66)
  q = ys(s1,s2,s3,s4,s5,s6,d12,d21,d23,d32,d13,d31,d44,d55,d66)

  plambdas, pInvariant = principalStress(p)   # no sort
  qlambdas, qInvariant = principalStress(q)   # no sort

  P1 = plambdas[0,:];  P2 = plambdas[1,:];  P3 = plambdas[2,:]
  Q1 = qlambdas[0,:];  Q2 = qlambdas[1,:];  Q3 = qlambdas[2,:]

  m2 = m/2.0;  m1 = 1.0/m;  m21 = m2-1.0
  P1Q1s = (P1-Q1)**2;  P1Q2s = (P1-Q2)**2;  P1Q3s = (P1-Q3)**2
  P2Q1s = (P2-Q1)**2;  P2Q2s = (P2-Q2)**2;  P2Q3s = (P2-Q3)**2
  P3Q1s = (P3-Q1)**2;  P3Q2s = (P3-Q2)**2;  P3Q3s = (P3-Q3)**2

  phi= P1Q1s**m2 + P1Q2s**m2 + P1Q3s**m2 + \
       P2Q1s**m2 + P2Q2s**m2 + P2Q3s**m2 + \
       P3Q1s**m2 + P3Q2s**m2 + P3Q3s**m2
  r  = (0.25*phi)**m1/sigma0 - 1.0

  if not Jac:
    return r.ravel()
  else:
    ln     = lambda x : np.log(x + 1.0e-32)

    drdphi =  (r+1.0)*m1/phi
    dphidm =( (P1Q1s**m2)*ln(P1Q1s) + (P1Q2s**m2)*ln(P1Q2s) + (P1Q3s**m2)*ln(P1Q3s) +
              (P2Q1s**m2)*ln(P2Q1s) + (P2Q2s**m2)*ln(P2Q2s) + (P2Q3s**m2)*ln(P2Q3s) +
              (P3Q1s**m2)*ln(P3Q1s) + (P3Q2s**m2)*ln(P3Q2s) + (P3Q3s**m2)*ln(P3Q3s) )*0.5
    jm     = drdphi*dphidm + (r+1.0)*ln(0.25*phi)*(-m1*m1)

    dPdc = principalStrs_Der(p, pInvariant, s1,s2,s3,s4,s5,s6)
    dQdd = principalStrs_Der(q, qInvariant, s1,s2,s3,s4,s5,s6)
    dphidP = m*np.array([ P1Q1s**m21*(P1-Q1) + P1Q2s**m21*(P1-Q2) + P1Q3s**m21*(P1-Q3),
                          P2Q1s**m21*(P2-Q1) + P2Q2s**m21*(P2-Q2) + P2Q3s**m21*(P2-Q3),
                          P3Q1s**m21*(P3-Q1) + P3Q2s**m21*(P3-Q2) + P3Q3s**m21*(P3-Q3) ])
    dphidQ = m*np.array([ P1Q1s**m21*(Q1-P1) + P2Q1s**m21*(Q1-P2) + P3Q1s**m21*(Q1-P3),
                          P1Q2s**m21*(Q2-P1) + P2Q2s**m21*(Q2-P2) + P3Q2s**m21*(Q2-P3),
                          P1Q3s**m21*(Q3-P1) + P2Q3s**m21*(Q3-P2) + P3Q3s**m21*(Q3-P3)])
    
    jc = drdphi*(dphidP[0]*dPdc[0]+dphidP[1]*dPdc[1]+dphidP[2]*dPdc[2])
    jd = drdphi*(dphidQ[0]*dQdd[0]+dphidQ[1]*dQdd[1]+dphidQ[2]*dQdd[2])
    return np.vstack((jc,jd, jm)).T

def KarafillisBoyceBasis(sigma0, c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,
                         b1, b2, a, alpha , sigmas, Jac=False):
  s1 = sigmas[0]; s2 = sigmas[1]; s3 = sigmas[2]
  s4 = sigmas[3]; s5 = sigmas[4]; s6 = sigmas[5]

  ks = lambda s1,s2,s3,s4,s5,s6,c11,c12,c13,c14,c15,c16: np.array( [
       ((c12+c13)*s1-c13*s2-c12*s3)/3.0,  ((c13+c11)*s2-c13*s1-c11*s3)/3.0,
       ((c11+c12)*s3-c12*s1-c11*s2)/3.0,  c14*s4, c15*s5, c16*s6    ])
  p = ks(s1,s2,s3,s4,s5,s6,c11,c12,c13,c14,c15,c16)
  q = ks(s1,s2,s3,s4,s5,s6,c21,c22,c23,c24,c25,c26)

  plambdas, pInvariant = principalStress(p)
  qlambdas, qInvariant = principalStress(q)

  P1 = plambdas[0,:];  P2 = plambdas[1,:];  P3 = plambdas[2,:]
  Q1 = qlambdas[0,:];  Q2 = qlambdas[1,:];  Q3 = qlambdas[2,:]

  b1h = b1/2.0; b1h1 = b1h-1.0; b2h = b2/2.0; b2h1 = b2h-1.0
  b1i = 1.0/b1; b2i = 1.0/b2
  ai = 1.0/a
  P2P3s = (P2-P3)**2;  Q1s = Q1**2
  P3P1s = (P3-P1)**2;  Q2s = Q2**2
  P1P2s = (P1-P2)**2;  Q3s = Q3**2

  phi10 = P2P3s**b1h + P3P1s**b1h + P1P2s**b1h
  phi20 = Q1s**b2h+Q2s**b2h+Q3s**b2h;  rb2 = 3.0**b2/(2.0**b2+2.0)
  phi1  = (0.5*phi10)**b1i
  phi2  = (rb2*phi20)**b2i

  Stress = alpha*phi1**a + (1.0-alpha)*phi2**a; EqStress = Stress**ai
  r  = EqStress/sigma0 - 1.0

  if not Jac:
    return r.ravel()
  else:
    ln     = lambda x : np.log(x + 1.0e-32)

    drds    = (r+1.0)*ai/Stress
    drdphi1 =  drds*     alpha *a*phi1**(a-1.0)
    drdphi2 =  drds*(1.0-alpha)*a*phi2**(a-1.0)
    dsda = alpha*phi1**a*ln(phi1) + (1.0-alpha)*phi2**a*ln(phi2)

    dphi1dphi10 = phi1/phi10/b1;  dphi2dphi20 = phi2/phi20/b2
    dphi1dP = np.array([ dphi1dphi10*b1*( P3P1s**b1h1*(P1-P3) + P1P2s**b1h1*(P1-P2)),
                         dphi1dphi10*b1*( P2P3s**b1h1*(P2-P3) + P1P2s**b1h1*(P2-P1)),
                         dphi1dphi10*b1*( P3P1s**b1h1*(P3-P1) + P2P3s**b1h1*(P3-P2)) ])
    dphi2dQ = np.array([ dphi2dphi20*b2*Q1s*b2h1*Q1,
                         dphi2dphi20*b2*Q2s*b2h1*Q2,
                         dphi2dphi20*b2*Q3s*b2h1*Q3 ])

    dPdc= principalStrs_Der(p, pInvariant, s1,s2,s3,s4,s5,s6, Karafillis=True)
    dQdc= principalStrs_Der(q, qInvariant, s1,s2,s3,s4,s5,s6, Karafillis=True)

    dphi10db1 = ( (P2P3s**b1h)*ln(P2P3s)+(P3P1s**b1h)*ln(P3P1s)+(P1P2s**b1h)*ln(P1P2s) )*0.5
    dphi20db2 = ( (P2P3s**b1h)*ln(P2P3s)+(P3P1s**b1h)*ln(P3P1s)+(P1P2s**b1h)*ln(P1P2s) )*0.5
    drb2db2   = rb2*ln(3.0) - rb2*ln(2.0)/(1.0+2.0**(1.0-b2))
    dphi1db1  =  phi1*ln(phi10)*(-b1i*b1i) + b1i*phi1/(0.5*phi10)* 0.5*dphi10db1
    dphi2db2  =  phi2*ln(phi20)*(-b2i*b2i) + b2i*phi2/(rb2*phi20)*(rb2*dphi20db2 + drb2db2*phi20)
    ja  = drds*dsda - (r+1.0)*ln(Stress)/a/a  #drda
    jb1 = drds * (     alpha *a*phi1**(a-1)) * dphi1db1
    jb2 = drds * ((1.0-alpha)*a*phi2**(a-1)) * dphi2db2
    jalpha = drds * (phi1**a - phi2**a)
    
    jc1 = drdphi1*(dphi1dP[0]*dPdc[0]+dphi1dP[1]*dPdc[1]+dphi1dP[2]*dPdc[2])
    jc2 = drdphi2*(dphi2dQ[0]*dQdc[0]+dphi2dQ[1]*dQdc[1]+dphi2dQ[2]*dQdc[2])

    return np.vstack((jc1,jc2,jb1,jb2,ja,jalpha)).T


fittingCriteria = {
  'tresca'         :{'func' : Tresca,
                     'num'  : 1,
                     'name' : 'Tresca',
                     'paras': 'Initial yield stress:',
                     'text' : '\nCoefficient of Tresca criterion:\nsigma0: ',
                     'error': 'The standard deviation error is: '
                    },
  'vonmises'       :{'func' : vonMises,
                     'num'  : 1,
                     'name' : 'Huber-Mises-Hencky(von Mises)',
                     'paras': 'Initial yield stress:',
                     'text' : '\nCoefficient of Huber-Mises-Hencky criterion:\nsigma0: ',
                     'error': 'The standard deviation error is: '
                    },
  'hosfordiso'     :{'func' : Hosford,
                     'num'  : 2,
                     'name' : 'Gerenal isotropic Hosford',
                     'paras': 'Initial yield stress, a:',
                     'text' : '\nCoefficients of Hosford criterion:\nsigma0, a: ',
                     'error': 'The standard deviation errors are: '
                    },
  'hosfordaniso'   :{'func' : generalHosford,
                     'num'  : 5,
                     'name' : 'Gerenal isotropic Hosford',
                     'paras': 'Initial yield stress, F, G, H, a:',
                     'text' : '\nCoefficients of Hosford criterion:\nsigma0, F, G, H, a: ',
                     'error': 'The standard deviation errors are: '
                    },
  'hill1948'       :{'func' : Hill1948,
                     'num'  : 6,
                     'name' : 'Hill1948',
                     'paras': 'Normalized [F, G, H, L, M, N]:',
                     'text' : '\nCoefficients of Hill1948 criterion:\n[F, G, H, L, M, N]:'+' '*16,
                     'error': 'The standard deviation errors are: '
                    },
  'hill1979'       :{'func' : Hill1979,
                     'num'  : 7,
                     'name' : 'Hill1979',
                     'paras': 'f,g,h,a,b,c,m:',
                     'text' : '\nCoefficients of Hill1979 criterion:\n f,g,h,a,b,c,m:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'drucker'        :{'func' : Drucker,
                     'num'  : 2,
                     'name' : 'Drucker',
                     'paras': 'Initial yield stress, C_D:',
                     'text' : '\nCoefficients of Drucker criterion:\nsigma0, C_D: ',
                     'error': 'The standard deviation errors are: '
                    },
  'gdrucker'       :{'func' : generalDrucker,
                     'num'  : 3,
                     'name' : 'General Drucker',
                     'paras': 'Initial yield stress, C_D, p:',
                     'text' : '\nCoefficients of Drucker criterion:\nsigma0, C_D, p: ',
                     'error': 'The standard deviation errors are: '
                    },
  'barlat1991iso'  :{'func' : Barlat1991iso,
                     'num'  : 2,
                     'name' : 'Barlat1991iso',
                     'paras': 'Initial yield stress, m:',
                     'text' : '\nCoefficients of isotropic Barlat 1991 criterion:\nsigma0, m:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'barlat1991aniso':{'func' : Barlat1991aniso,
                     'num'  : 8,
                     'name' : 'Barlat1991aniso',
                     'paras': 'Initial yield stress, a, b, c, f, g, h, m:',
                     'text' : '\nCoefficients of anisotropic Barlat 1991 criterion:\nsigma0, a, b, c, f, g, h, m:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'bbc2003'        :{'func' : BBC2003,
                     'num'  : 9,
                     'name' : 'Banabic-Balan-Comsa 2003',
                     'paras': 'Initial yield stress, a, b, c, d, e, f, g, k:',
                     'text' : '\nCoefficients of anisotropic Barlat 1991 criterion:\nsigma0, a, b, c, d, e, f, g, k:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'bbc2005'        :{'func' : BBC2005,
                     'num'  : 9,'err':np.inf,
                     'name' : 'Banabic-Balan-Comsa 2003',
                     'paras': 'a, b, L ,M, N, P, Q, R, k:',
                     'text' : '\nCoefficients of Banabic-Balan-Comsa 2005 criterion: a, b, L ,M, N, P, Q, R, k:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'Cazacu_Barlat2D':{'func' : Cazacu_Barlat2D,
                     'num'  : 11,
                     'name' : 'Cazacu Barlat for plain stress',
                     'paras': 'a1,a2,a3,a6; b1,b2,b3,b4,b5,b10; c:',
                     'text' : '\nCoefficients of Cazacu Barlat yield criterion for plane stress: \
                               \n a1,a2,a3,a6; b1,b2,b3,b4,b5,b10; c:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'Cazacu_Barlat3D':{'func' : Cazacu_Barlat3D,
                     'num'  : 18,
                     'name' : 'Cazacu Barlat',
                     'paras': 'a1,a2,a3,a4,a5,a6; b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11; c:',
                     'text' : '\nCoefficients of Cazacu Barlat yield criterion for plane stress: \
                               \n a1,a2,a3,a4,a5,a6; b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11; c\n',
                     'error': 'The standard deviation errors are: '
                    },
  'yld200418p'     :{'func' : Yld200418p,
                     'num'  : 20,
                     'name' : 'Yld200418p',
                     'paras': 'Equivalent stress,c12,c21,c23,c32,c31,c13,c44,c55,c66,d12,d21,d23,d32,d31,d13,d44,d55,d66,m:',
                     'text' : '\nCoefficients of Yld2004-18p yield criterion: \
                               \n Y, c12,c21,c23,c32,c31,c13,c44,c55,c66,d12,d21,d23,d32,d31,d13,d44,d55,d66,m\n',
                     'error': 'The standard deviation errors are: '
                    },
  'karafillis'     :{'func' : KarafillisBoyce,
                     'num'  : 16,
                     'name' : 'Yld200418p',
                     'paras': 'c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,b1,b2,a,alpha',
                     'text' : '\nCoefficients of Karafillis-Boyce yield criterion: \
                               \n c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,b1,b2,a,alpha\n',
                     'error': 'The standard deviation errors are: '
                    }
                  }

for key in fittingCriteria.keys():
  if 'num' in fittingCriteria[key].keys(): 
    fittingCriteria[key]['bound']=[(None,None)]*fittingCriteria[key]['num']
    fittingCriteria[key]['guess']=np.ones(fittingCriteria[key]['num'],'d')

thresholdParameter = ['totalshear','equivalentStrain']

#---------------------------------------------------------------------------------------------------
class Loadcase():
#---------------------------------------------------------------------------------------------------
  '''
     Class for generating load cases for the spectral solver
  '''

# ------------------------------------------------------------------
  def __init__(self,finalStrain,incs,time,ND=3,RD=1,nSet=1,dimension=3,vegter=False):
    print('using the random load case generator')
    self.finalStrain = finalStrain
    self.incs = incs
    self.time = time
    self.ND = ND
    self.RD = RD
    self.nSet = nSet
    self.dimension = dimension
    self.vegter = vegter
    self.NgeneratedLoadCases = 0
    if self.vegter:
      self.vegterLoadcase = self._vegterLoadcase()

  def getLoadcase(self,number):
    if self.dimension == 3:
      print 'generate random 3D load case'
      return self._getLoadcase3D()
    else:
      if self.vegter is True:
        print 'generate load case for Vegter'
        return self._getLoadcase2dVegter(number)
      else:
        print 'generate random 2D load case'
        return self._getLoadcase2dRandom()

  def getLoadcase3D(self):
    self.NgeneratedLoadCases+=1
    defgrad=['*']*9
    stress =[0]*9
    values=(np.random.random_sample(9)-.5)*self.finalStrain*2

    main=np.array([0,4,8])
    np.random.shuffle(main)
    for i in main[:2]:                                                                              # fill 2 out of 3 main entries
      defgrad[i]=1.+values[i]
      stress[i]='*'
    for off in [[1,3,0],[2,6,0],[5,7,0]]:                                                           # fill 3 off-diagonal pairs of defgrad (1 or 2 entries)
      off=np.array(off)
      np.random.shuffle(off)
      for i in off[0:2]:
        if i != 0: 
          defgrad[i]=values[i]
          stress[i]='*'

    return 'f '+' '.join(str(c) for c in defgrad)+\
          ' p '+' '.join(str(c) for c in stress)+\
          ' incs %s'%self.incs+\
          ' time %s'%self.time

  def _getLoadcase2dVegter(self,number): #for a 2D simulation, I would use this generator before switching to a random 2D generator
    NDzero=[[1,2,3,6],[1,3,5,7],[2,5,6,7]] # no deformation / * for stress
    # biaxial f1 = f2
    # shear   f1 = -f2
    # unixaial f1 , f2 =0
    # plane strain f1 , s2 =0
    # modulo to get one out of 4
    stress =['*', '*', '0']*3
    defgrad = self.vegterLoadcase[number-1]

    return 'f '+' '.join(str(c) for c in defgrad)+\
          ' p '+' '.join(str(c) for c in stress)+\
          ' incs %s'%self.incs+\
          ' time %s'%self.time

  def _vegterLoadcase(self):
    '''
      generate the stress points for Vegter criteria
    '''
    theta = np.linspace(0.0,np.pi/2.0,self.nSet)
    f = [0.0, 0.0, '*']*3;  loadcase = []
    for i in xrange(self.nSet*4): loadcase.append(f)

    # more to do for F
    F = np.array([  [[1.1, 0.1], [0.1, 1.1]],              # uniaxial tension
                    [[1.1, 0.1], [0.1, 1.1]],              # shear
                    [[1.1, 0.1], [0.1, 1.1]],              # eq-biaxial
                    [[1.1, 0.1], [0.1, 1.1]],              # eq-biaxial
                 ])
    for i,t in enumerate(theta):
      R = np.array([np.cos(t), np.sin(t), -np.sin(t), np.cos(t)]).reshape(2,2)
      for j in xrange(4): 
        loadcase[i*4+j][0],loadcase[i*4+j][1],loadcase[i*4+j][3],loadcase[i*4+j][4] = np.dot(R.T,np.dot(F[j],R)).reshape(4)
    return loadcase

  def _getLoadcase2dRandom(self):
    '''
      generate random stress points for 2D tests
    '''
    self.NgeneratedLoadCases+=1
    defgrad=['0', '0', '*']*3
    stress =['*', '*', '0']*3
    defgrad[0],defgrad[1],defgrad[3],defgrad[4] = (np.random.random_sample(4)-.5)*self.finalStrain*2.0 + np.eye(2).reshape(4)

    return 'f '+' '.join(str(c) for c in defgrad)+\
          ' p '+' '.join(str(c) for c in stress)+\
          ' incs %s'%self.incs+\
          ' time %s'%self.time
  def _defgradScale(self, defgrad, finalStrain):
    '''
    '''
    defgrad0 = (np.array([ 0.0 if i is '*' else i for i in defgrad ]))
    det0 = 1.0 - numpy.linalg.det(defgrad0.reshape(3,3))
    if defgrad0[0] == 0.0: defgrad0[0] = det0/(defgrad0[4]*defgrad0[8]-defgrad0[5]*defgrad0[7])
    if defgrad0[4] == 0.0: defgrad0[4] = det0/(defgrad0[0]*defgrad0[8]-defgrad0[2]*defgrad0[6])
    if defgrad0[8] == 0.0: defgrad0[8] = det0/(defgrad0[0]*defgrad0[4]-defgrad0[1]*defgrad0[3])
    strain   = np.dot(defgrad0.reshape(3,3).T,defgrad0.reshape(3,3)) - np.eye(3)
    eqstrain = 2.0/3.0*np.sqrt( 1.5*(strain[0][0]**2+strain[1][1]**2+strain[2][2]**2) + 
                                3.0*(strain[0][1]**2+strain[1][2]**2+strain[2][0]**2) )
    r = finalStrain*1.25/eqstrain
#    if r>1.0: defgrad =( np.array([i*r if i is not '*' else i for i in defgrad]))


#---------------------------------------------------------------------------------------------------
class Criterion(object):
#---------------------------------------------------------------------------------------------------
  '''
     Fitting to certain criterion
  '''
  def __init__(self,name='worst'):
    self.name = name
    self.results = fittingCriteria

    if self.name.lower() not in map(str.lower, self.results.keys()):
      raise Exception('no suitable fitting criterion selected')
    else:
      print('fitting to the %s criterion'%name)

  def fit(self,stress):
    global fitResults

    nameCriterion = self.name.lower()
    criteriaClass = fittingCriteria[nameCriterion]['func']
    numParas      = fittingCriteria[nameCriterion]['num']
    textParas     = fittingCriteria[nameCriterion]['text'] + formatOutput(numParas)
    textError     = fittingCriteria[nameCriterion]['error']+ formatOutput(numParas,'%-14.8f')+'\n'
    bounds        = fittingCriteria[nameCriterion]['bound']  # Default bounds, no bound
    guess0        = fittingCriteria[nameCriterion]['guess']  # Default initial guess, depends on bounds
    criteria      = criteriaClass(0.0)
    if fitResults == [] : initialguess = guess0
    else                : initialguess = np.array(fitResults[-1])

    weight = get_weight(np.shape(stress)[1])
    ydata  = np.zeros(np.shape(stress)[1])
    try:
      popt, pcov, infodict, errmsg, ierr = \
         leastsqBound (criteria.fun,  initialguess,        args=(ydata,stress),
                       bounds=bounds, Dfun=criteria.jac,   full_output=True)
      if ierr not in [1, 2, 3, 4]:
        raise RuntimeError("Optimal parameters not found: " + errmsg)
      if (len(ydata) > len(initialguess)) and pcov is not None:
        s_sq = (criteria.fun(popt, *(ydata,stress))**2).sum()/(len(ydata)-len(initialguess))
        pcov = pcov * s_sq
      perr = np.sqrt(np.diag(pcov))
      fitResults.append(popt.tolist())

      print (textParas%array2tuple(popt))
      print (textError%array2tuple(perr))
      print('Number of function calls =', infodict['nfev'])
    except Exception as detail:
      print detail
      pass


#---------------------------------------------------------------------------------------------------
class myThread (threading.Thread):
#---------------------------------------------------------------------------------------------------
  '''
     Runner class
  '''
  def __init__(self, threadID):
    threading.Thread.__init__(self)
    self.threadID = threadID
  def run(self):
    s.acquire()
    conv=converged()
    s.release()
    while not conv:
      doSim(4.,self.name)
      s.acquire()
      conv=converged()
      s.release()

def doSim(delay,thread):
  
  s.acquire()
  me=loadcaseNo()
  if not os.path.isfile('%s.load'%me):
    print('generating loadcase for sim %s from %s'%(me,thread))
    f=open('%s.load'%me,'w')
    f.write(myLoad.getLoadcase(me))
    f.close()
    s.release()
  else: s.release()
  
  s.acquire()
  if not os.path.isfile('%s_%i.spectralOut'%(options.geometry,me)):
    print('starting simulation %s from %s'%(me,thread))
    s.release()
    execute('DAMASK_spectral -g %s -l %i'%(options.geometry,me))
  else: s.release()
   
  s.acquire()
  if not os.path.isfile('./postProc/%s_%i.txt'%(options.geometry,me)):
    print('starting post processing for sim %i from %s'%(me,thread))
    s.release()
    try:
      execute('postResults --cr f,p --co totalshear %s_%i.spectralOut'%(options.geometry,me))
    except:
      execute('postResults --cr f,p %s_%i.spectralOut'%(options.geometry,me))
    execute('addCauchy ./postProc/%s_%i.txt'%(options.geometry,me))
    execute('addStrainTensors -l -v ./postProc/%s_%i.txt'%(options.geometry,me))
    execute('addMises -s Cauchy -e ln(V) ./postProc/%s_%i.txt'%(options.geometry,me))
  else: s.release()

  s.acquire()
  print('-'*10)
  print('reading values for sim %i from %s'%(me,thread))
  s.release()

  refFile = open('./postProc/%s_%i.txt'%(options.geometry,me))
  table = damask.ASCIItable(refFile)
  table.head_read()
  if options.fitting =='equivalentStrain':
    thresholdKey = 'Mises(ln(V))'
  elif options.fitting =='totalshear':
    thresholdKey = 'totalshear'
  s.acquire()
  for l in [thresholdKey,'1_Cauchy']:
    if l not in table.labels: print '%s not found'%l
  s.release()
  table.data_readArray(['%i_Cauchy'%(i+1) for i in xrange(9)]+[thresholdKey]+['%i_ln(V)'%(i+1) for i in xrange(9)])

  line = 0
  lines = np.shape(table.data)[0]
  yieldStress     = np.empty((int(options.yieldValue[2]),6),'d')
  deformationRate = np.empty((int(options.yieldValue[2]),6),'d')
  for i,threshold in enumerate(np.linspace(options.yieldValue[0],options.yieldValue[1],options.yieldValue[2])):
    while line < lines:
      if table.data[line,9]>= threshold:
        upper,lower = table.data[line,9],table.data[line-1,9]                                       # values for linear interpolation
        stress = np.array(table.data[line-1,0:9] * (upper-threshold)/(upper-lower) + \
                          table.data[line  ,0:9] * (threshold-lower)/(upper-lower)).reshape(3,3)    # linear interpolation of stress values
        dstrain= np.array(table.data[line,10:] - table.data[line-1,10:]).reshape(3,3)

        yieldStress[i,0]= stress[0,0]; yieldStress[i,1]=stress[1,1]; yieldStress[i,2]=stress[2,2]
        yieldStress[i,3]=(stress[0,1] + stress[1,0])/2.0     #   0  3  5
        yieldStress[i,4]=(stress[1,2] + stress[2,1])/2.0     #   *  1  4  yieldStress
        yieldStress[i,5]=(stress[2,0] + stress[0,2])/2.0     #   *  *  2

#       D*dt = 0.5(L+L^T)*dt = 0.5*d(lnF + lnF^T) = dlnV
        deformationRate[i,0]= dstrain[0,0]; deformationRate[i,1]=dstrain[1,1]; deformationRate[i,2]=dstrain[2,2]
        deformationRate[i,3]=(dstrain[0,1] + dstrain[1,0])/2.0     #   0  3  5
        deformationRate[i,4]=(dstrain[1,2] + dstrain[2,1])/2.0     #   *  1  4
        deformationRate[i,5]=(dstrain[2,0] + dstrain[0,2])/2.0     #   *  *  2
        break
      else:
        line+=1

  s.acquire()
  global stressAll, strainAll
  print('number of yield points of sim %i: %i'%(me,len(yieldStress)))
  print('starting fitting for sim %i from %s'%(me,thread))
  try:
    for i in xrange(int(options.yieldValue[2])):
      stressAll[i]=np.append(stressAll[i], yieldStress[i]/unitGPa)
      strainAll[i]=np.append(strainAll[i], deformationRate[i])
      myFit.fit(stressAll[i].reshape(len(stressAll[i])//6,6).transpose())
  except Exception as detail:
    print('could not fit for sim %i from %s'%(me,thread))
    print detail
    s.release()
    return
  s.release()

def loadcaseNo():
  global N_simulations
  N_simulations+=1
  return N_simulations

def converged():
  global N_simulations
  if N_simulations < options.max:
    return False
  else:
    return True

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Performs calculations with various loads on given geometry file and fits yield surface.

""", version=string.replace(scriptID,'\n','\\n')
)
# maybe make an option to specifiy if 2D/3D fitting should be done?
parser.add_option('-l','--load' ,    dest='load', type='float', nargs=3,
                                     help='load: final strain; increments; time %default', metavar='float int float')
parser.add_option('-g','--geometry', dest='geometry', type='string',
                                     help='name of the geometry file [%default]', metavar='string')
parser.add_option('-c','--criterion',     dest='criterion', choices=fittingCriteria.keys(),
                                     help='criterion for stopping simulations [%default]', metavar='string')
parser.add_option('-f','--fitting',       dest='fitting', choices=thresholdParameter,
                                     help='yield criterion [%default]', metavar='string')
parser.add_option('-y','--yieldvalue',    dest='yieldValue', type='float', nargs=3,
                                     help='yield points: start; end; count %default', metavar='float float int')
parser.add_option('--min',           dest='min', type='int',
                                     help='minimum number of simulations [%default]', metavar='int')
parser.add_option('--max',           dest='max', type='int',
                                     help='maximum number of iterations [%default]',  metavar='int')
parser.add_option('-t','--threads',       dest='threads', type='int',
                                     help='number of parallel executions [%default]',  metavar='int')
parser.add_option('-d','--dimension',     dest='dimension', type='int',
                                     help='dimension of the virtual test [%default]',  metavar='int')
parser.add_option('-v', '--vegter',  dest='vegter', action='store_true',
                  help='Vegter criteria [%default]')                                     
parser.set_defaults(min        = 12)
parser.set_defaults(max        = 30)
parser.set_defaults(threads    = 4)
parser.set_defaults(yieldValue = (0.002,0.004,2))
parser.set_defaults(load       = (0.010,100,100.0))
parser.set_defaults(criterion  = 'worst')
parser.set_defaults(fitting    = 'totalshear')
parser.set_defaults(geometry   = '20grains16x16x16')
parser.set_defaults(dimension  = 3)
parser.set_defaults(vegter     = 'False')


options = parser.parse_args()[0]

if not os.path.isfile(options.geometry+'.geom'):
  parser.error('geometry file %s.geom not found'%options.geometry)
if not os.path.isfile('material.config'):
  parser.error('material.config file not found')
if options.threads<1:
  parser.error('invalid number of threads %i'%options.threads)
if options.min<0:
  parser.error('invalid minimum number of simulations %i'%options.min)
if options.max<options.min:
  parser.error('invalid maximum number of simulations (below minimum)')
if options.yieldValue[0]>options.yieldValue[1]:
  parser.error('invalid yield start (below yield end)')
if options.yieldValue[2] != int(options.yieldValue[2]):
  parser.error('count must be an integer')
if not os.path.isfile('numerics.config'):
  print('numerics.config file not found')

if not os.path.isfile('material.config'):
  print('material.config file not found')

if options.vegter is True: 
  options.dimension = 2
unitGPa = 10.e8
N_simulations=0
fitResults = []
s=threading.Semaphore(1)

stressAll=[np.zeros(0,'d').reshape(0,0) for i in xrange(int(options.yieldValue[2]))]
strainAll=[np.zeros(0,'d').reshape(0,0) for i in xrange(int(options.yieldValue[2]))]
myLoad = Loadcase(options.load[0],options.load[1],options.load[2],
                  nSet = 10, dimension = options.dimension, vegter = options.vegter)
myFit = Criterion(options.criterion)

threads=[]

for i in range(options.threads):
  threads.append(myThread(i))
  threads[i].start()

for i in range(options.threads):
  threads[i].join()

print 'finished fitting to yield criteria'
