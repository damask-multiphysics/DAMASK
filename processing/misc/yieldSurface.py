#!/usr/bin/python
# -*- coding: UTF-8 no BOM -*-

import threading,time,os,subprocess,string,sys
import numpy as np
from optparse import OptionParser
import damask
from damask.util import leastsqBound

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

def runFit(exponent, eqStress, dimension, criterion):
  global s, threads, myFit, myLoad
  global fitResults, fitErrors, fitResidual, stressAll, strainAll
  global N_simulations, Guess, dDim
  
  fitResults = []; fitErrors = []; fitResidual = []; Guess = []; threads=[]
  dDim = dimension - 3
  Nparas = len(fitCriteria[criterion]['bound'][dDim])
  nExpo = fitCriteria[criterion]['nExpo']

  if exponent > 0.0:                                                                                # User defined exponents
    Nparas = Nparas-nExpo
    fitCriteria[criterion]['bound'][dDim] = fitCriteria[criterion]['bound'][dDim][:Nparas]

  for i in xrange(Nparas):
    temp = fitCriteria[criterion]['bound'][dDim][i]
    if fitCriteria[criterion]['bound'][dDim][i] == (None,None):
      Guess.append(1.0)
    else: 
      g = (temp[0]+temp[1])/2.0
      if g == 0: g = temp[1]*0.5
      Guess.append(g)

  N_simulations=0
  s=threading.Semaphore(1)
  myLoad = Loadcase(options.load[0],options.load[1],options.load[2],
                    nSet = 10, dimension = dimension, vegter = options.vegter)

  stressAll= [np.zeros(0,'d').reshape(0,0) for i in xrange(int(options.yieldValue[2]))]
  strainAll= [np.zeros(0,'d').reshape(0,0) for i in xrange(int(options.yieldValue[2]))]

  myFit = Criterion(exponent,eqStress, dimension, criterion)
  for t in range(options.threads):
    threads.append(myThread(t))
    threads[t].start()

  for t in range(options.threads):
    threads[t].join()
  damask.util.croak('Residuals')
  damask.util.croak(fitResidual)

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

def principalStress(p):
  I = invariant(p)

  I1s3I2= (I[0]**2 - 3.0*I[1])**0.5
  numer = 2.0*I[0]**3 - 9.0*I[0]*I[1] + 27.0*I[2]
  denom = 2.0*I1s3I2**3
  cs    = numer/denom

  phi   = np.arccos(cs)/3.0
  t1    = I[0]/3.0; t2 = 2.0/3.0*I1s3I2
  return np.array( [t1 + t2*np.cos(phi), 
                    t1 + t2*np.cos(phi+np.pi*2.0/3.0),
                    t1 + t2*np.cos(phi+np.pi*4.0/3.0)])

def principalStrs_Der(p, (s1, s2, s3, s4, s5, s6), dim, Karafillis=False):
  '''
    Derivative of principal stress with respect to stress
  '''
  third  = 1.0/3.0
  third2 = 2.0*third

  I = invariant(p)
  I1s3I2= np.sqrt(I[0]**2 - 3.0*I[1])
  numer = 2.0*I1**3 - 9.0*I[0]*I[1] + 27.0*I[2]
  denom = 2.0*I1s3I2**3
  cs    = numer/denom
  phi   = np.arccos(cs)/3.0

  dphidcs   = -third/np.sqrt(1.0 - cs**2)
  dcsddenom = 0.5*numer*(-1.5)*I1s3I2**(-5.0)
  dcsdI1    = (6.0*I1**2 - 9.0*I2)*denom + dcsddenom*(2.0*I1)
  dcsdI2    = (          - 9.0*I1)*denom + dcsddenom*(-3.0)
  dcsdI3    = 27.0*denom
  dphidI1, dphidI2, dphidI3 = dphidcs*dcsdI1, dphidcs*dcsdI2, dphidcs*dcsdI3

  dI1s3I2dI1 =   I1/I1s3I2
  dI1s3I2dI2 = -1.5/I1s3I2
  tcoeff  = third2*I1s3I2
  
  dSidIj = lambda theta : ( tcoeff*(-np.sin(theta))*dphidI1 + third2*dI1s3I2dI1*np.cos(theta) + third, 
                            tcoeff*(-np.sin(theta))*dphidI2 + third2*dI1s3I2dI2*np.cos(theta), 
                            tcoeff*(-np.sin(theta))*dphidI3)
  dSdI = np.array([dSidIj(phi),dSidIj(phi+np.pi*2.0/3.0),dSidIj(phi+np.pi*4.0/3.0)]) # i=1,2,3; j=1,2,3
  
# calculate the derivation of principal stress with regards to the anisotropic coefficients  
  one    = np.ones_like(s1); zero = np.zeros_like(s1); num  = len(s1)
  dIdp = np.array([[one,    one,    one,    zero,   zero,   zero], 
                   [p[1]+p[2], p[2]+p[0], p[0]+p[1], -2.0*p[3], -2.0*p[4], -2.0*p[5]],
                   [p[1]*p[2]-p[4]**2, p[2]*p[0]-p[5]**2, p[0]*p[1]-p[3]**2, 
                    -2.0*p[3]*p[2]+2.0*p[4]*p[5], -2.0*p[4]*p[0]+2.0*p[5]*p[3], -2.0*p[5]*p[1]+2.0*p[3]*p[4]] ]) 
  if Karafillis:
    dpdc = np.array([[zero,s1-s3,s1-s2], [s2-s3,zero,s2-s1], [s3-s2,s3-s1,zero]])/3.0
    dSdp = np.array([np.dot(dSdI[:,:,i],dIdp[:,:,i]).T for i in xrange(num)]).T
    if dim == 2:
      temp = np.vstack([dSdp[:,3]*s4]).T.reshape(num,1,3).T
    else: 
      temp = np.vstack([dSdp[:,3]*s4,dSdp[:,4]*s5,dSdp[:,5]*s6]).T.reshape(num,3,3).T

    return np.concatenate((np.array([np.dot(dSdp[:,0:3,i], dpdc[:,:,i]).T for i in xrange(num)]).T,
                           temp), axis=1)
  else:
    if dim == 2:
      dIdc=np.array([[-dIdp[i,0]*s2, -dIdp[i,1]*s1, -dIdp[i,1]*s3,
                      -dIdp[i,2]*s2, -dIdp[i,2]*s1, -dIdp[i,0]*s3,
                       dIdp[i,3]*s4 ] for i in xrange(3)])
    else:
      dIdc=np.array([[-dIdp[i,0]*s2, -dIdp[i,1]*s1, -dIdp[i,1]*s3,
                      -dIdp[i,2]*s2, -dIdp[i,2]*s1, -dIdp[i,0]*s3,
                       dIdp[i,3]*s4,  dIdp[i,4]*s5,  dIdp[i,5]*s6 ] for i in xrange(3)])
    return np.array([np.dot(dSdI[:,:,i],dIdc[:,:,i]).T for i in xrange(num)]).T

def invariant(sigmas):
    I=np.zeros(3)
    s11,s22,s33,s12,s23,s31 = sigmas
    I[0] = s11 + s22 + s33
    I[1] = s11*s22 + s22*s33 + s33*s11 - s12**2 - s23**2 - s31**2
    I[2]  = s11*s22*s33 + 2.0*s12*s23*s31 - s12**2*s33 - s23**2*s11 - s31**2*s22
    return I

def math_ln(x):
  return np.log(x + 1.0e-32)

def sym6to33(sigma6):
  ''' Shape the symmetric stress tensor(6,1) into (3,3) '''
  sigma33 = np.empty((3,3))
  sigma33[0,0] = sigma6[0]; sigma33[1,1] = sigma6[1]; sigma33[2,2] = sigma6[2];
  sigma33[0,1] = sigma6[3]; sigma33[1,0] = sigma6[3]
  sigma33[1,2] = sigma6[4]; sigma33[2,1] = sigma6[4]
  sigma33[2,0] = sigma6[5]; sigma33[0,2] = sigma6[5]
  return sigma33


class Criteria(object):
  '''
    residuum of anisotropic Barlat 1991 yield criterion (eq. 2.37)
  '''
  def __init__(self, criterion, uniaxialStress,exponent, dimension):
    self.stress0 = uniaxialStress
    if exponent < 0.0:                                                                              # Fitting exponent m
      self.mFix = [False, exponent]
    else:                                                                                           # fixed exponent m
      self.mFix = [True, exponent]
    self.func     = fitCriteria[criterion]['func']
    self.criteria = criterion
    self.dim      = dimension
  def fun(self, paras, ydata, sigmas):
    return self.func(self.stress0, paras, sigmas,self.mFix,self.criteria,self.dim)
  def jac(self, paras, ydata, sigmas):
    return self.func(self.stress0, paras, sigmas,self.mFix,self.criteria,self.dim,Jac=True)

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
    damask.util.croak('Warning: the number of stress points is uncorrect, stress points of %s are missing in set %i'%(
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

def Tresca(eqStress, paras, sigmas, mFix, criteria, dim, Jac = False):
  '''
    Tresca yield criterion
    the fitted parameters is: paras(sigma0)
    eqStress, mFix, criteria, dim are invalid input
  '''
  if not Jac:
    lambdas = principalStresses(sigmas)
    r = np.amax(np.array([abs(lambdas[2,:]-lambdas[1,:]),\
                          abs(lambdas[1,:]-lambdas[0,:]),\
                          abs(lambdas[0,:]-lambdas[2,:])]),0) - paras
    return r.ravel()
  else:
    return -np.ones(len(sigmas))

def Cazacu_Barlat(eqStress, paras, sigmas, mFix, criteria, dim, Jac = False):
  '''
    Cazacu-Barlat (CB) yield criterion
    the fitted parameters are:
      a1,a2,a3,a6; b1,b2,b3,b4,b5,b10; c for plane stress
      a1,a2,a3,a4,a5,a6; b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11; c: for general case
    mFix are invalid input
  '''
  s11,s22,s33,s12,s23,s31 = sigmas
  if dim == 2:
    (a1,a2,a3,a4), (b1,b2,b3,b4,b5,b10), c = paras[0:4],paras[4:10],paras[10]
    a5 = a6 = b6 = b7 = b8 = b9 = b11 = 0.0
    s33 = s23 = s31 = np.zeros_like(s11)
  else:
    (a1,a2,a3,a4,a5,a6), (b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11), c = paras[0:6],paras[6:17],paras[17]

  s1_2, s2_2, s3_2, s12_2, s23_2, s31_2 = np.array([s11,s22,s33,s12,s23,s31])**2
  s1_3, s2_3, s3_3, s123, s321 = s11*s1_2, s22*s2_2, s33*s3_2,s11*s22*s33, s12*s23*s31
  d12_2,d23_2,d31_2 = (s11-s22)**2, (s22-s33)**2, (s33-s11)**2
  
  J20 = ( a1*d12_2 + a2*d23_2 + a3*d31_2 )/6.0 + a4*s12_2 + a5*s23_2 + a6*s31_2
  J30 = ( (b1    +b2    )*s1_3  + (b3    +b4    )*s2_3  + ( b1+b4-b2      +  b1+b4-b3     )*s3_3 )/27.0- \
        ( (b1*s22+b2*s33)*s1_2  + (b3*s33+b4*s11)*s2_2  + ((b1+b4-b2)*s11 + (b1+b4-b3)*s22)*s3_2 )/9.0 + \
        ( (b1+b4)*s123/9.0 + b11*s321 )*2.0 - \
        ( ( 2.0*b9 *s22 - b8*s33  - (2.0*b9 -b8)*s11 )*s31_2 + 
          ( 2.0*b10*s33 - b5*s22  - (2.0*b10-b5)*s11 )*s12_2 +
          ( (b6+b7)*s11 - b6*s22  - b7*s33           )*s23_2
        )/3.0
  f0  = J20**3 - c*J30**2
  r   = f0**(1.0/6.0)*np.sqrt(3.0)/eqStress

  if not Jac:
    return (r - 1.0).ravel()
  else:
    drdf = r/f0/6.0
    dj2, dj3 = drdf*3.0*J20**2, -drdf*2.0*J30*c
    jc   = -drdf*J30**2

    ja1,ja2,ja3 =  dj2*d12_2/6.0, dj2*d23_2/6.0, dj2*d31_2/6.0
    ja4,ja5,ja6 =  dj2*s12_2,     dj2*s23_2,     dj2*s31_2
    jb1 = dj3*( (s1_3 + 2.0*s3_3)/27.0 - s22*s1_2/9.0 - (s11+s22)*s3_2/9.0 + s123/4.5 )
    jb2 = dj3*( (s1_3 -     s3_3)/27.0 - s33*s1_2/9.0 +  s11     *s3_2/9.0 )
    jb3 = dj3*( (s2_3 -     s3_3)/27.0 - s33*s2_2/9.0 +  s22     *s3_2/9.0 )
    jb4 = dj3*( (s2_3 + 2.0*s3_3)/27.0 - s11*s2_2/9.0 - (s11+s22)*s3_2/9.0 + s123/4.5 )

    jb5, jb10 = dj3*(s22 - s11)*s12_2/3.0,  dj3*(s11 - s33)*s12_2/1.5
    jb6, jb7  = dj3*(s22 - s11)*s23_2/3.0,  dj3*(s33 - s11)*s23_2/3.0
    jb8, jb9  = dj3*(s33 - s11)*s31_2/3.0,  dj3*(s11 - s22)*s31_2/1.5
    jb11      = dj3*s321*2.0
    if dim == 2:
      return np.vstack((ja1,ja2,ja3,ja4,jb1,jb2,jb3,jb4,jb5,jb10,jc)).T
    else:
      return np.vstack((ja1,ja2,ja3,ja4,ja5,ja6,jb1,jb2,jb3,jb4,jb5,jb6,jb7,jb8,jb9,jb10,jb11,jc)).T

def Drucker(eqStress, paras, sigmas, mFix, criteria, dim, Jac = False):
  '''
    Drucker yield criterion
    the fitted parameters are 
      sigma0, C_D for Drucker(p=1);
      sigma0, C_D, p for general Drucker
    eqStress, mFix are invalid inputs
  '''
  if criteria == 'drucker':
    sigma0, C_D= paras
    p = 1.0
  else:
    sigma0, C_D = paras[0:2]
    if mFix[0]: p = mFix[1]
    else:       p = paras[-1]
  I = invariant(sigmas)
  J = np.zeros([3]) 
  J[1]  = I[0]**2/3.0  - I[1]
  J[2]  = I[0]**3/13.5 - I[0]*I[1]/3.0 + I[2]
  J2_3p = J[1]**(3.0*p)      
  J3_2p = J[2]**(2.0*p)
  left  = J2_3p - C_D*J3_2p
  r     = left**(1.0/(6.0*p))*3.0**0.5/sigma0

  if not Jac:
    return (r - 1.0).ravel()
  else:
    drdl = r/left/(6.0*p)
    if criteria == 'drucker':
      return np.vstack((-r/sigma0, -drdl*J3_2p)).T
    else:
      dldp = 3.0*J2_3p*math_ln(J2) - 2.0*C_D*J3_2p*math_ln(J3)
      jp   = drdl*dldp + r*math_ln(left)/(-6.0*p*p)
      
      if mFix[0]: return np.vstack((-r/sigma0, -drdl*J3_2p)).T
      else:       return np.vstack((-r/sigma0, -drdl*J3_2p, jp)).T

def Hill1948(eqStress, paras, sigmas, mFix, criteria, dim, Jac = False):
  '''
    Hill 1948 yield criterion
    the fitted parameters are:
      F, G, H, L, M, N for 3D
      F, G, H, N       for 2D
    eqStress, mFix, criteria are invalid input
  '''
  s11,s22,s33,s12,s23,s31 = sigmas
  if dim == 2:      # plane stress
    jac = np.array([ s22**2, s11**2, (s11-s22)**2, 2.0*s12**2])
  else:             # general case
    jac = np.array([(s22-s33)**2,(s33-s11)**2,(s11-s22)**2, 2.0*s23**2,2.0*s31**2,2.0*s12**2])
  if not Jac:
    return (np.dot(paras,jac)/2.0-0.5).ravel()
  else:
    return jac.T

def Hill1979(eqStress,paras, sigmas, mFix, criteria, dim, Jac = False):
  '''
    Hill 1979 yield criterion
    the fitted parameters are: f,g,h,a,b,c,m
    criteria are invalid input
  '''
  if mFix[0]: m = mFix[1]
  else:       m = paras[-1]

  coeff  = paras[0:6]
  s1,s2,s3 = principalStresses(sigmas)
  # s= principalStresses(sigmas)
  diffs  = np.array([s2-s3, s3-s1, s1-s2, 2.0*s1-s2-s3, 2.0*s2-s3-s1, 2.0*s3-s1-s2])**2 
  #diffs  = np.array([s[1]-s[2], s[2]-s[0], etc ...  s1-s2, 2.0*s1-s2-s3, 2.0*s2-s3-s1, 2.0*s3-s1-s2])**2 
  diffsm = diffs**(m/2.0)
  left   = np.dot(coeff,diffsm)
  r = (0.5*left)**(1.0/m)/eqStress #left = base**mi

  if not Jac:
    return (r-1.0).ravel()
  else:
    drdl, dldm = r/left/m, np.dot(coeff,diffsm*math_ln(diffs))*0.5
    jm = drdl*dldm + r*math_ln(0.5*left)*(-1.0/m/m) #/(-m**2) 

    if mFix[0]: return np.vstack((drdl*diffsm)).T
    else:       return np.vstack((drdl*diffsm, jm)).T

def Hosford(eqStress, paras, sigmas, mFix, criteria, dim, Jac = False):
  '''
    Hosford family criteria
    the fitted parameters are:
      von Mises: sigma0
      Hershey: (1) sigma0, a, when a is not fixed; (2) sigma0, when a is fixed
      general Hosford: (1) F,G,H, a, when a is not fixed; (2) F,G,H, when a is fixed
  '''

  if criteria == 'vonmises':
    sigma0 = paras
    coeff  = np.ones(3)
    a = 2.0
  elif criteria == 'hershey':
    sigma0 = paras[0]
    coeff  = np.ones(3)
    if mFix[0]: a = mFix[1]
    else:       a = paras[1]
  else:
    sigma0 = eqStress
    coeff  = paras[0:3]
    if mFix[0]: a = mFix[1]
    else:       a = paras[3]

  s1,s2,s3 = principalStresses(sigmas)
  diffs  = np.array([s2-s3, s3-s1, s1-s2])**2
  diffsm = diffs**(a/2.0)
  left   = np.dot(coeff,diffsm)
  r      = (0.5*left)**(1.0/a)/sigma0

  if not Jac:
    return (r-1.0).ravel()
  else:
    if criteria == 'vonmises': # von Mises
      return -r/sigma0
    else:
      drdl, dlda = r/left/a, np.dot(coeff,diffsm*math_ln(diffs))*0.5
      ja = drdl*dlda + r*math_ln(0.5*left)*(-1.0/a/a)
      if criteria == 'hershey':  # Hershey
        if mFix[0]: return -r/sigma0
        else:       return np.vstack((-r/sigma0, ja)).T
      else:                                        # Anisotropic Hosford
        if mFix[0]: return np.vstack((drdl*diffsm)).T
        else:       return np.vstack((drdl*diffsm, ja)).T

def Barlat1989(eqStress, paras, sigmas, mFix, criteria, dim, Jac=False):
  '''
    Barlat-Lian 1989 yield criteria
    the fitted parameters are:
      Anisotropic: a, h, p, m; m is optional
  ''' 
  a, h, p = paras[0:3]
  if mFix[0]: m = mFix[1]
  else:       m = paras[-1]

  c = 2.0-a
  s11,s22,s12 = sigmas[0], sigmas[1], sigmas[3]
  k1,k2 = 0.5*(s11 + h*s22), (0.25*(s11 - h*s22)**2 + (p*s12)**2)**0.5
  fs  = np.array([ (k1+k2)**2, (k1-k2)**2, 4.0*k2**2 ]); fm = fs**(m/2.0)                
  left = np.dot(np.array([a,a,c]),fm)
  r = (0.5*left)**(1.0/m)/eqStress

  if not Jac:
    return (r-1.0).ravel()
  else:
    dk1dh = 0.5*s22
    dk2dh, dk2dp = 0.25*(s11-h*s22)*(-s22)/k2, p*s12**2/k2
    dlda,  dldc  = fm[0]+fm[1], fm[2]
    fm1 = fs**(m/2.0-1.0)*m
    dldk1, dldk2 = a*fm1[0]*(k1+k2)+a*fm1[1]*(k1-k2), a*fm1[0]*(k1+k2)-a*fm1[1]*(k1-k2)+c*fm1[2]*k2*4.0
    drdl,  drdm  = r/m/left, r*math_ln(0.5*left)*(-1.0/m/m) 
    dldm = np.dot(np.array([a,a,c]),fm*math_ln(fs))*0.5

    ja = drdl*dlda 
    jh,jp = drdl*(dldk1*dk1dh + dldk2*dk2dh), drdl*dldk2*dk2dp
    jm    = drdl*dldm + drdm

    if mFix[0]: return np.vstack((ja,jc,jh,jp)).T
    else:       return np.vstack((ja,jc,jh,jp,jm)).T

def Barlat1991(eqStress, paras, sigmas, mFix, criteria, dim, Jac=False):
  '''
    Barlat 1991 criteria
    the fitted parameters are:
      Anisotropic: a, b, c, f, g, h, m for 3D
                   a, b, c, h, m for plane stress
    m is optional
  '''
  if dim == 2: coeff = paras[0:4]    # plane stress  
  else:        coeff = paras[0:6]    # general case
  if mFix[0]:  m = mFix[1]
  else:        m = paras[-1]

  s11,s22,s33,s12,s23,s31 = sigmas
  if dim == 2:
    dXdx = np.array([s22,-s11,s11-s22,s12])
    A,B,C,H = np.array(coeff)[:,None]*dXdx; F=G=0.0
  else:
    dXdx = np.array([s22-s33,s33-s11,s11-s22,s23,s31,s12])
    A,B,C,F,G,H = np.array(coeff)[:,None]*dXdx

  I2 = (F*F + G*G + H*H)/3.0+ ((A-C)**2+(C-B)**2+(B-A)**2)/54.0
  I3 = (C-B)*(A-C)*(B-A)/54.0 + F*G*H - ((C-B)*F*F + (A-C)*G*G + (B-A)*H*H)/6.0
  phi1 = np.arccos(I3/I2**1.5)/3.0 + np.pi/6.0; absc1 = 2.0*np.abs(np.cos(phi1))
  phi2 = phi1                      + np.pi/3.0; absc2 = 2.0*np.abs(np.cos(phi2))
  phi3 = phi2                      + np.pi/3.0; absc3 = 2.0*np.abs(np.cos(phi3))
  left = ( absc1**m + absc2**m + absc3**m )
  r    = (0.5*left)**(1.0/m)*np.sqrt(3.0*I2)/eqStress

  if not Jac:
    return (r - 1.0).ravel()
  else:
    dfdl = r/left/m
    jm = r*math_ln(0.5*left)*(-1.0/m/m) + dfdl*0.5*(
         absc1**m*math_ln(absc1) + absc2**m*math_ln(absc2) + absc3**m*math_ln(absc3) ) 

    da,db,dc = (2.0*A-B-C)/18.0, (2.0*B-C-A)/18.0, (2.0*C-A-B)/18.0
    if dim == 2:
      dI2dx = np.array([da, db, dc, H])/1.5*dXdx
      dI3dx = np.array([ da*(B-C) + (H**2-G**2)/2.0, 
                         db*(C-A) + (F**2-H**2)/2.0, 
                         dc*(A-B) + (G**2-F**2)/2.0,
                         (G*F + (A-B))*H ])/3.0*dXdx
    else:
      dI2dx = np.array([da, db, dc, F,G,H])/1.5*dXdx
      dI3dx = np.array([ da*(B-C) + (H**2-G**2)/2.0, 
                         db*(C-A) + (F**2-H**2)/2.0, 
                         dc*(A-B) + (G**2-F**2)/2.0,
                        (H*G*3.0 + (B-C))*F, 
                        (F*H*3.0 + (C-A))*G, 
                        (G*F*3.0 + (A-B))*H ])/3.0*dXdx
    darccos = -1.0/np.sqrt(1.0 - I3**2/I2**3)

    dfdcos = lambda phi : dfdl*m*(2.0*abs(np.cos(phi)))**(m-1.0)*np.sign(np.cos(phi))*(-np.sin(phi)/1.5)

    dfdthe= (dfdcos(phi1) + dfdcos(phi2) + dfdcos(phi3)) 
    dfdI2, dfdI3 = dfdthe*darccos*I3*(-1.5)*I2**(-2.5)+r/2.0/I2, dfdthe*darccos*I2**(-1.5)

    if mFix[0]: return np.vstack((dfdI2*dI2dx + dfdI3*dI3dx)).T
    else:       return np.vstack((dfdI2*dI2dx + dfdI3*dI3dx, jm)).T

def BBC2000(eqStress, paras, sigmas, mFix, criteria, dim, Jac=False):
  '''
    BBC2000 yield criterion
    the fitted parameters are 
    d,e,f,g, b,c,a, k;  k is optional
    criteria are invalid input
  '''
  d,e,f,g, b,c,a= paras[0:7]
  if mFix[0]: k = mFix[1]
  else:       k = paras[-1]

  s11,s22,s12 = sigmas[0], sigmas[1], sigmas[3]
  k2  = 2.0*k;  k1 = k - 1.0
  M,N,P,Q,R = d+e, e+f, (d-e)/2.0, (e-f)/2.0, g**2
  Gamma =    M*s11 + N*s22
  Psi   = ( (P*s11 + Q*s22)**2 + s12**2*R )**0.5

  l1, l2, l3  = b*Gamma + c*Psi, b*Gamma - c*Psi, 2.0*c*Psi
  l1s,l2s,l3s = l1**2,           l2**2,           l3**2
                       
  left = a*l1s**k + a*l2s**k + (1-a)*l3s**k
  r = left**(1.0/k2)/eqStress
  if not Jac:
    return (r - 1.0).ravel()
  else:
    drdl,drdk = r/left/k2, r*math_ln(left)*(-1.0/k2/k)
    dldl1,dldl2,dldl3 = a*k2*(l1s**k1)*l1, a*k2*(l2s**k1)*l2, (1-a)*k2*(l3s**k1)*l3
    dldGama, dldPsi = (dldl1 + dldl2)*b, (dldl1 - dldl2 + 2.0*dldl3)*c
    temp = (P*s11 + Q*s22)/Psi
    dPsidP, dPsidQ, dPsidR = temp*s11, temp*s22, 0.5*s12**2/Psi
    dlda = l1s**k + l2s**k - l3s**k
    dldb = dldl1*Gamma + dldl2*Gamma
    dldc = dldl1*Psi   - dldl2*Psi + dldl3*2.0*Psi
    dldk = a*math_ln(l1s)*l1s**k + a*math_ln(l2s)*l2s**k + (1-a)*math_ln(l3s)*l3s**k

    J = drdl*np.array([dldGama*s11+dldPsi*dPsidP*0.5, dldGama*(s11+s22)+dldPsi*(-dPsidP+dPsidQ)*0.5, #jd,je
                       dldGama*s22-dldPsi*dPsidQ*0.5, dldPsi*dPsidR*2.0*g,                           #jf,jg
                       dldb, dldc, dlda])                                                            #jb,jc,ja
    if mFix[0]: return np.vstack(J).T
    else:       return np.vstack((J, drdl*dldk + drdk)).T


def BBC2003(eqStress, paras, sigmas, mFix, criteria, dim, Jac=False):
  '''
    BBC2003 yield criterion
    the fitted parameters are 
    M,N,P,Q,R,S,T,a, k;  k is optional
    criteria are invalid input
  '''
  M,N,P,Q,R,S,T,a = paras[0:8]
  if mFix[0]: k = mFix[1]
  else:       k = paras[-1]

  s11,s22,s12 = sigmas[0], sigmas[1], sigmas[3]
  k2  = 2.0*k;  k1 = k - 1.0
  Gamma  =   0.5 *  (s11 + M*s22)
  Psi    = ( 0.25*(N*s11 - P*s22)**2 + Q*Q*s12**2 )**0.5
  Lambda = ( 0.25*(R*s11 - S*s22)**2 + T*T*s12**2 )**0.5

  l1, l2, l3  = Gamma + Psi, Gamma - Psi, 2.0*Lambda
  l1s,l2s,l3s = l1**2,        l2**2,        l3**2
  left = a*l1s**k + a*l2s**k + (1-a)*l3s**k
  r = left**(1.0/k2)/eqStress
  if not Jac:
    return (r - 1.0).ravel()
  else:
    drdl,drdk = r/left/k2, r*math_ln(left)*(-1.0/k2/k)
    dldl1,dldl2,dldl3 = a*k2*(l1s**k1)*l1, a*k2*(l2s**k1)*l2, (1-a)*k2*(l3s**k1)*l3

    dldGamma, dldPsi, dldLambda  = dldl1+dldl2, dldl1-dldl2, 2.0*dldl3
    temp = 0.25/Psi*(N*s11 - P*s22)
    dPsidN, dPsidP, dPsidQ = s11*temp, -s22*temp, Q*s12**2/Psi
    temp = 0.25/Lambda*(R*s11 - S*s22)
    dLambdadR, dLambdadS, dLambdadT = s11*temp, -s22*temp, T*s12**2/Psi
    dldk = a*math_ln(l1s)*l1s**k + a*math_ln(l2s)*l2s**k + (1-a)*math_ln(l3s)*l3s**k

    J = drdl * np.array([dldGamma*s22*0.5,                                              #jM
                         dldPsi*dPsidN,       dldPsi*dPsidP,       dldPsi*dPsidQ,       #jN, jP, jQ
                         dldLambda*dLambdadR, dldLambda*dLambdadS, dldLambda*dLambdadT, #jR, jS, jT
                         l1s**k + l2s**k - l3s**k ])                                    #ja

    if mFix[0]: return np.vstack(J).T
    else :      return np.vstack((J, drdl*dldk+drdk)).T

def BBC2005(eqStress, paras, sigmas, mFix, criteria, dim, Jac=False):
  '''
    BBC2005 yield criterion
    the fitted parameters are 
    a, b, L ,M, N, P, Q, R, k;  k is optional
    criteria are invalid input
  '''
  a,b,L, M, N, P, Q, R = paras[0:8]
  if mFix[0]: k = mFix[1]
  else:       k = paras[-1]

  s11 = sigmas[0]; s22 = sigmas[1]; s12 = sigmas[3]
  k2  = 2.0*k
  Gamma  =    L*s11 + M*s22
  Lambda = ( (N*s11 - P*s22)**2 + s12**2 )**0.5
  Psi    = ( (Q*s11 - R*s22)**2 + s12**2 )**0.5

  l1  = Lambda + Gamma; l2  = Lambda - Gamma; l3  = Lambda + Psi; l4  = Lambda - Psi
  l1s = l1**2;          l2s = l2**2;          l3s = l3**2;        l4s = l4**2
  left = a*l1s**k + a*l2s**k + b*l3s**k + b*l4s**k
  sBar = left**(1.0/k2);  r = sBar/eqStress - 1.0
  if not Jac:
    return r.ravel()
  else:
    ln   = lambda x : np.log(x + 1.0e-32)
    expo = 0.5/k;  k1 = k-1.0

    dsBardl = expo*sBar/left/eqStress
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

    J = dsBardl * np.array( [
      l1s**k+l2s**k, l3s**k+l4s**k,dldGama*s11,dldGama*s22,dldLambda*dLambdadN,
      dldLambda*dLambdadP, dldPsi*dPsidQ, dldPsi*dPsidR])

    if mFix[0]: return np.vstack(J).T
    else :      return np.vstack(J, dldk+dsBarde*dedk).T

def Yld2000(eqStress, paras, sigmas, mFix, criteria, dim, Jac=False):
  '''
    C: c11,c22,c66  c12=c21=1.0 JAC NOT PASS
    D: d11,d12,d21,d22,d66 
  '''
  C,D = paras[0:3], paras[3:8]
  if mFix[0]: m = mFix[1]
  else:       m = paras[-1]

  s11, s22, s12 = sigmas[0],sigmas[1],sigmas[3]
  X = np.array([ 2.0*C[0]*s11-C[0]*s22, 2.0*C[1]*s22-C[1]*s11, 3.0*C[2]*s12 ])/3.0 # a1,a2,a7
  Y = np.array([ (8.0*D[2]-2.0*D[0]-2.0*D[3]+2.0*D[1])*s11 + (4.0*D[3]-4.0*D[1]-4.0*D[2]+    D[0])*s22,
                 (4.0*D[0]-4.0*D[2]-4.0*D[1]+    D[3])*s11 + (8.0*D[1]-2.0*D[3]-2.0*D[0]+2.0*D[2])*s22,
                  9.0*D[4]*s12 ])/9.0

  def priStrs((sx,sy,sxy)):
    temp = np.sqrt( (sx-sy)**2 + 4.0*sxy**2 )
    return 0.5*(sx+sy + temp), 0.5*(sx+sy - temp)
  m2 = m/2.0; m21 = m2 - 1.0
  (X1,X2), (Y1,Y2) = priStrs(X), priStrs(Y) # Principal values of X, Y
  phi1s, phi21s, phi22s  = (X1-X2)**2, (2.0*Y2+Y1)**2, (2.0*Y1+Y2)**2
  phi1,  phi21,  phi22 = phi1s**m2, phi21s**m2, phi22s**m2
  left = phi1 + phi21 + phi22
  r = (0.5*left)**(1.0/m)/eqStress

  if not Jac:
    return (r-1.0).ravel()
  else:
    drdl,  drdm  = r/m/left, r*math_ln(0.5*left)*(-1.0/m/m) #/(-m*m) 
    dldm = ( phi1*math_ln(phi1s) + phi21*math_ln(phi21s) + phi22*math_ln(phi22s) )*0.5
    zero = np.zeros_like(s11); num  = len(s11)
    def dPrincipalds((X1,X2,X12)): 
      # the derivative of principla with regards to stress
      temp   = 1.0/np.sqrt( (X1-X2)**2 + 4.0*X12**2 )
      dP1dsi = 0.5*np.array([ 1.0+temp*(X1-X2), 1.0-temp*(X1-X2),  temp*4.0*X12])
      dP2dsi = 0.5*np.array([ 1.0-temp*(X1-X2), 1.0+temp*(X1-X2), -temp*4.0*X12])
      return np.array([dP1dsi, dP2dsi])
    
    dXdXi, dYdYi = dPrincipalds(X), dPrincipalds(Y)
    dXidC  = np.array([ [ 2.0*s11-s22,        zero,    zero ],  #dX11dC
                        [        zero, 2.0*s22-s11,    zero ],  #dX22dC
                        [        zero,        zero, 3.0*s12 ] ])/3.0  #dX12dC
    dYidD  = np.array([ [ -2.0*s11+    s22,  2.0*s11-4.0*s22,  8.0*s11-4.0*s22, -2.0*s11+4.0*s22,    zero ],  #dY11dD
                        [  4.0*s11-2.0*s22, -4.0*s11+8.0*s22, -4.0*s11+2.0*s22,      s11-2.0*s22,    zero ],  #dY22dD
                        [             zero,             zero,             zero,             zero, 9.0*s12 ] ])/9.0  #dY12dD
    
    dXdC=np.array([np.dot(dXdXi[:,:,i], dXidC[:,:,i]).T for i in xrange(num)]).T
    dYdD=np.array([np.dot(dYdYi[:,:,i], dYidD[:,:,i]).T for i in xrange(num)]).T
    
    dldX = m*np.array([ phi1s**m21*(X1-X2),           phi1s**m21*(X2-X1)])
    dldY = m*np.array([phi21s**m21*(2.0*Y2+Y1) + 2.0*phi22s**m21*(2.0*Y1+Y2), \
                       phi22s**m21*(2.0*Y1+Y2) + 2.0*phi21s**m21*(2.0*Y2+Y1) ])
    jC = drdl*np.array([np.dot(dldX[:,i], dXdC[:,:,i]) for i in xrange(num)]).T
    jD = drdl*np.array([np.dot(dldY[:,i], dYdD[:,:,i]) for i in xrange(num)]).T

    jm   = drdl*dldm + drdm
    if mFix[0]: return np.vstack((jC,jD)).T
    else:       return np.vstack((jC,jD,jm)).T

def Yld200418p(eqStress, paras, sigmas, mFix, criteria, dim, Jac=False):
  '''
    Yld2004-18p yield criterion
    the fitted parameters are 
      C: c12,c21,c23,c32,c31,c13,c44,c55,c66; D: d12,d21,d23,d32,d31,d13,d44,d55,d66 for 3D
      C: c12,c21,c23,c32,c31,c13,c44; D: d12,d21,d23,d32,d31,d13,d44 for 2D
    and m, m is optional
    criteria are invalid input
  '''
  if dim == 2: C,D = np.append(paras[0:7],[0.0,0.0]), np.append(paras[7:14],[0.0,0.0])
  else:        C,D = paras[0:9], paras[9:18]
  if mFix[0]: m = mFix[1]
  else:       m = paras[-1]

  sv = (sigmas[0] + sigmas[1] + sigmas[2])/3.0
  sdev = np.vstack((sigmas[0:3]-sv,sigmas[3:6]))
  ys = lambda sdev, C: np.array([-C[0]*sdev[1]-C[5]*sdev[2], -C[1]*sdev[0]-C[2]*sdev[2], 
                                 -C[4]*sdev[0]-C[3]*sdev[1],  C[6]*sdev[3], C[7]*sdev[4], C[8]*sdev[5]])
  p,q = ys(sdev, C), ys(sdev, D)
  pLambdas, qLambdas = principalStress(p), principalStress(q)   # no sort

  m2 = m/2.0; x3 = xrange(3); num = len(sv)
  PiQj  = np.array([(pLambdas[i,:]-qLambdas[j,:]) for i in x3 for j in x3])
  QiPj  = np.array([(qLambdas[i,:]-pLambdas[j,:]) for i in x3 for j in x3]).reshape(3,3,num)
  PiQjs = PiQj**2
  left  = np.sum(PiQjs**m2,axis=0)
  r = (0.25*left)**(1.0/m)/eqStress

  if not Jac:
    return (r - 1.0).ravel()
  else:
    drdl, drdm = r/m/left, r*math_ln(0.25*left)*(-1.0/m/m)
    dldm = np.sum(PiQjs**m2*math_ln(PiQjs),axis=0)*0.5
    dPdc, dQdd = principalStrs_Der(p, sdev, dim), principalStrs_Der(q, sdev, dim)
    PiQjs3d = ( PiQjs**(m2-1.0) ).reshape(3,3,num)
    dldP = -m*np.array([np.diag(np.dot(PiQjs3d[:,:,i], QiPj   [:,:,i])) for i in xrange(num)]).T
    dldQ =  m*np.array([np.diag(np.dot(QiPj   [:,:,i], PiQjs3d[:,:,i])) for i in xrange(num)]).T

    jm = drdl*dldm + drdm
    jc = drdl*np.sum([dldP[i]*dPdc[i] for i in x3],axis=0)
    jd = drdl*np.sum([dldQ[i]*dQdd[i] for i in x3],axis=0)

    if mFix[0]: return np.vstack((jc,jd)).T
    else:       return np.vstack((jc,jd,jm)).T

def KarafillisBoyce(eqStress, paras, sigmas, mFix, criteria, dim, Jac=False):
  '''
    Karafillis-Boyce
    the fitted parameters are 
      c11,c12,c13,c14,c15,c16,c,m for 3D
      c11,c12,c13,c14,c,m for plane stress
    0<c<1, m are optional
    criteria are invalid input
  '''
  ks = lambda (s1,s2,s3,s4,s5,s6),(c1,c2,c3,c4,c5,c6): np.array( [
         ((c2+c3)*s1-c3*s2-c2*s3)/3.0,  ((c3+c1)*s2-c3*s1-c1*s3)/3.0,
         ((c1+c2)*s3-c2*s1-c1*s2)/3.0,  c4*s4, c5*s5, c6*s6 ])
  if dim == 2: C1,c = np.append(paras[0:4],[0.0,0.0]), paras[4]
  else:        C1,c = paras[0:6], paras[6]
  if mFix[0]:  m = mFix[1]
  else:        m = paras[-1]  # Karafillis-Boyce

  p= ks(sigmas, C1)
  plambdas = principalStress(p)
  reci_m, m2, rm2, m1 = 1.0/m, m/2.0, 3.0**m/(2.0**(m-1.0)+1.0), m-1.0

  difP  = np.array([ plambdas[0]-plambdas[1], plambdas[1]-plambdas[2], plambdas[2]-plambdas[0] ])
  difPs = difP**2; difPm1 = difPs**(m2-1.0)
  Ps    = plambdas**2

  phi1, phi2 = np.sum(difPs**m2, axis = 0), np.sum(Ps**m2, axis = 0)
  left = (1.0-c)*phi1+ c*rm2*phi2
  r    = (0.5*left)**reci_m/eqStress

  if not Jac:
    return (r-1.0).ravel()
  else:
    drdl, drdm = r*reci_m/left,  -r*math_ln(0.5*left)*reci_m*reci_m
    dldm = (1.0-c)*np.sum(difPs**m2*math_ln(difPs), axis=0)*0.5 + \
            rm2*c *np.sum(   Ps**m2*math_ln(Ps),    axis=0)*0.5 + \
            rm2*c *phi2* ( np.log(3.0) - 2.0**m1/(2.0**m1 + 1.0)*np.log(2.0) )
    dphi1dP = m*np.array([ difPm1[0]*difP[0] - difPm1[2]*difP[2], 
                           difPm1[1]*difP[1] - difPm1[0]*difP[0], 
                           difPm1[2]*difP[2] - difPm1[1]*difP[1] ])
    dphi2dP = m*plambdas*Ps**(m2-1.0)

    dPdc = principalStrs_Der(p, sigmas, dim, Karafillis=True)
    dldP = (1.0-c)*dphi1dP + c*dphi2dP*rm2

    jm  = drdl * dldm + drdm  #drda*(-1.0/m/m)
    jc1 = drdl * np.sum([dldP[i]*dPdc[i] for i in xrange(3)],axis=0)
    jc  = drdl * (-phi1 + rm2*phi2)

    if mFix[0]: return np.vstack((jc1,jc)).T
    else:       return np.vstack((jc1,jc,jm)).T



fitCriteria = {
  'tresca'         :{'name':  'Tresca',
                     'func':   Tresca,
                     'nExpo':  0,'err':np.inf,
                     'dimen':  3,
                     'bound':  [[(None,None)]],
                     'labels': [['sigma0']],
                    },
  'vonmises'       :{'name':  'Huber-Mises-Hencky',
                     'func' :  Hosford,
                     'nExpo':  0,'err':np.inf,
                     'dimen':  3,
                     'bound':  [[(None,None)]],
                     'labels': [['sigma0']],
                    },
  'hershey'        :{'name':   'Hershey',
                     'func':   Hosford,
                     'nExpo':  1,'err':np.inf,
                     'dimen':  3, 
                     'bound':  [[(None,None)]+[(1.0,8.0)]],
                     'labels': [['sigma0','a']],
                    },
  'hosford'        :{'name':   'General Hosford',
                     'func':   Hosford,
                     'nExpo':  1,'err':np.inf,
                     'dimen':  3,
                     'bound':  [[(0.0,2.0)]*3+[(1.0,8.0)] ],
                     'labels': [['F','G','H','a']],
                    },
  'hill1948'       :{'name':   'Hill 1948',
                     'func':   Hill1948,
                     'nExpo':  0,'err':np.inf,
                     'dimen':  3,
                     'bound':  [[(None,None)]*6, [(None,None)]*4 ],
                     'labels': [['F','G','H','L','M','N'],['F','G','H','N']],
                    },
  'hill1979'       :{'name':   'Hill 1979',
                     'func':   Hill1979,
                     'nExpo':  1,'err':np.inf,
                     'dimen':  3,
                     'bound':  [[(-2.0,2.0)]*6+[(1.0,8.0)] ],
                     'labels': [['f','g','h','a','b','c','m']],
                    },
  'drucker'        :{'name':   'Drucker',
                     'func':   Drucker,
                     'nExpo':  0,'err':np.inf,
                     'dimen':  3,
                     'bound':  [[(None,None)]+[(-3.375, 2.25)]],
                     'labels': [['sigma0','C_D']],
                    },
  'gdrucker'       :{'name':   'General Drucker',
                     'func':   Drucker,
                     'nExpo':  1,'err':np.inf,
                     'dimen':  3,
                     'bound':  [[(None,None)]+[(-3.375, 2.25)]+[(1.0,8.0)] ],
                     'labels': [['sigma0','C_D', 'p']],
                    },
  'barlat1989'     :{'name':   'Barlat 1989',
                     'func':   Barlat1989,
                     'nExpo':  1,'err':np.inf,
                     'dimen':  2,
                     'bound':  [[(-3.0,3.0)]*4+[(1.0,8.0)] ],
                     'labels': [['a','c','h','f', 'm']],
                    },
  'barlat1991'     :{'name':   'Barlat 1991',
                     'func':   Barlat1991,
                     'nExpo':  1,'err':np.inf,
                     'dimen':  3,
                     'bound':  [[(-2,2)]*6+[(1.0,8.0)], [(-2,2)]*4+[(1.0,8.0)]],
                     'labels': [['a','b','c','f','g','h','m'],['a','b','c','f','m']],
                    },
  'bbc2000'        :{'name':   'Banabic-Balan-Comsa 2000',
                     'func':   BBC2000,
                     'nExpo':  1,'err':np.inf,
                     'dimen':  2,
                     'bound':  [[(None,None)]*7+[(1.0,8.0)]], #[(None,None)]*6+[(0.0,1.0)]+[(1.0,9.0)],
                     'labels': [['d','e','f','g','b','c','a','k']],
                    },
  'bbc2003'        :{'name':   'Banabic-Balan-Comsa 2003',
                     'func' :  BBC2003,
                     'nExpo':  1,'err':np.inf,
                     'dimen':  2,
                     'bound':  [[(None,None)]*8+[(1.0,8.0)]], #[(None,None)]*7+[(0.0,1.0)]+[(1.0,9.0)],
                     'labels': [['M','N','P','Q','R','S','T','a','k']],
                    },
  'bbc2005'        :{'name':   'Banabic-Balan-Comsa 2005',
                     'func' :  BBC2005,
                     'nExpo':  1,'err':np.inf,
                     'dimen':  2,
                     'bound':  [[(None,None)]*8+[(1.0,8.0)] ], #[(None,None)]*6+[(0.0,1.0)]*2+[(1.0,9.0)],
                     'labels': [['L','M','N','P','Q','R','a','b','k']],
                    },
  'cazacu'         :{'name':   'Cazacu Barlat',
                     'func':   Cazacu_Barlat,
                     'nExpo':  0,'err':np.inf,
                     'dimen':  3,
                     'bound':  [[(None,None)]*16+[(-2.5,2.5)]+[(None,None)]],
                     'labels': [['a1','a2','a3','a4','a5','a6', 'b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','b11', 'c'],
                                ['a1','a2','a3','a6', 'b1','b2','b3','b4','b5','b10', 'c']],
                    },
  'yld2000'        :{'name':   'Yld2000-2D',
                     'func':   Yld2000,
                     'nExpo':  1,'err':np.inf,
                     'dimen':  2,
                     'bound':  [[(None,None)]*8+[(1.0,8.0)]],
                     'labels': [['a1','a2','a7','a3','a4','a5','a6','a8','m']],
                    },
  'yld200418p'     :{'name':   'Yld2004-18p',
                     'func' :  Yld200418p,
                     'nExpo':  1,'err':np.inf,
                     'dimen':  3,
                     'bound':  [[(None,None)]*18+[(1.0,8.0)], [(None,None)]*14+[(1.0,8.0)]],
                     'labels': [['c12','c21','c23','c32','c31','c13','c44','c55','c66','d12','d21','d23','d32','d31','d13','d44','d55','d66','m'],
                                ['c12','c21','c23','c32','c31','c13','c44','d12','d21','d23','d32','d31','d13','d44','m']],
                    },
  'karafillis'     :{'name':   'Karafillis-Boyce',
                     'func' :  KarafillisBoyce,
                     'nExpo':  1,'err':np.inf,
                     'dimen':  3,
                     'bound':  [[(None,None)]*6+[(0.0,1.0)]+[(1.0,8.0)], [(None,None)]*4+[(0.0,1.0)]+[(1.0,8.0)]],
                     'labels': [['c11','c12','c13','c14','c15','c16','c','m'],
                                ['c11','c12','c13','c14','c15','c16','c','m']],
                    },
                  }

thresholdParameter = ['totalshear','equivalentStrain']

#---------------------------------------------------------------------------------------------------
class Loadcase():
#---------------------------------------------------------------------------------------------------
  '''
     Class for generating load cases for the spectral solver
  '''

# ------------------------------------------------------------------
  def __init__(self,finalStrain,incs,time,ND=3,RD=1,nSet=1,dimension=3,vegter=False):
    damask.util.croak('Using the random load case generator')
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
      damask.util.croak('Generate random 3D load case')
      return self._getLoadcase3D()
    else:
      if self.vegter is True:
        damask.util.croak('Generate load case for Vegter')
        return self._getLoadcase2dVegter(number)
      else:
        damask.util.croak('Generate random 2D load case')
        return self._getLoadcase2dRandom()

  def _getLoadcase3D(self):
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
    ratio = self._defgradScale(defgrad)
    for i in [0,4,8]:
      if defgrad[i] != '*': defgrad[i] = (defgrad[i]-1.0)*ratio + 1.0
    for i in [1,2,3,5,6,7]:
      if defgrad[i] != '*': defgrad[i] = defgrad[i]*ratio

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
    # for i,t in enumerate(theta):
      # R = np.array([np.cos(t), np.sin(t), -np.sin(t), np.cos(t)]).reshape(2,2)
      # for j in xrange(4): 
        # loadcase[i*4+j][0],loadcase[i*4+j][1],loadcase[i*4+j][3],loadcase[i*4+j][4] = np.dot(R.T,np.dot(F[j],R)).reshape(4)
    # return loadcase

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
  def _defgradScale(self, defgrad):
    '''
    '''
    def fill_star(a,b):
      if   a != '*' and b != '*': return a,b
      elif a == '*' and b != '*': return b,b
      elif a != '*' and b == '*': return a,a
      else                      : return 0.0,0.0
    defgrad0 = defgrad[:]
    defgrad0[1],defgrad0[3] = fill_star(defgrad[1], defgrad[3])
    defgrad0[2],defgrad0[6] = fill_star(defgrad[2], defgrad[6])
    defgrad0[5],defgrad0[7] = fill_star(defgrad[5], defgrad[7])
    for i in [0,4,8]:
      if defgrad0[i] == '*': defgrad0[i] = 0.0
    det0 = 1.0 - np.linalg.det(np.array(defgrad0).reshape(3,3))
    if defgrad0[0] == 0.0: defgrad0[0] = det0/(defgrad0[4]*defgrad0[8]-defgrad0[5]*defgrad0[7])
    if defgrad0[4] == 0.0: defgrad0[4] = det0/(defgrad0[0]*defgrad0[8]-defgrad0[2]*defgrad0[6])
    if defgrad0[8] == 0.0: defgrad0[8] = det0/(defgrad0[0]*defgrad0[4]-defgrad0[1]*defgrad0[3])
    strain   = 0.5*(np.dot(np.array(defgrad0).reshape(3,3).T,np.array(defgrad0).reshape(3,3)) - np.eye(3)) #Green Strain
    eqstrain = 2.0/3.0*np.sqrt( 1.5*(strain[0][0]**2+strain[1][1]**2+strain[2][2]**2) + 
                                3.0*(strain[0][1]**2+strain[1][2]**2+strain[2][0]**2) )
    ratio =  self.finalStrain*1.05/eqstrain
    return max(ratio,1.0)

#---------------------------------------------------------------------------------------------------
class Criterion(object):
#---------------------------------------------------------------------------------------------------
  '''
     Fitting to certain criterion
  '''
  def __init__(self, exponent, uniaxial, dimension, label='vonmises'):
    self.name = label
    self.expo = exponent
    self.uniaxial= uniaxial
    self.dimen   = dimension
    self.results = fitCriteria

    if self.name.lower() not in map(str.lower, self.results.keys()):
      raise Exception('No suitable fitting criterion selected')
    else:
      damask.util.croak('Fitting to the %s criterion'%fitCriteria[self.name]['name'])
  
  def report_labels(self):
    if len(fitCriteria[self.name]['labels']) > 1 and self.dimen == 2:
      return fitCriteria[self.name]['labels'][1]
    else:
      return fitCriteria[self.name]['labels'][0]

  def report_name(self):
    return fitCriteria[self.name]['name']

  def fit(self,stress):
    global fitResults; fitErrors; fitResidual
    if options.exponent > 0.0: nExponent = nExpo
    else:                      nExponent = 0
    nameCriterion = self.name.lower()
    criteria      = Criteria(nameCriterion,self.uniaxial,self.expo, self.dimen)
    bounds        = fitCriteria[nameCriterion]['bound'][dDim]                                       # Default bounds, no bound
    guess0        = Guess                                                                           # Default initial guess, depends on bounds

    if fitResults == [] : initialguess = guess0
    else                : initialguess = np.array(fitResults[-1])

    ydata  = np.zeros(np.shape(stress)[1])
    try:
      popt, pcov, infodict, errmsg, ierr = \
         leastsqBound (criteria.fun,  initialguess,        args=(ydata,stress),
                       bounds=bounds, Dfun=criteria.jac,   full_output=True)
      if ierr not in [1, 2, 3, 4]:
        raise RuntimeError("Optimal parameters not found: "+errmsg)
      else:
        residual = criteria.fun(popt, ydata, stress)
        fitResidual.append(np.linalg.norm(residual)/np.sqrt(len(residual)))
      if (len(ydata) > len(initialguess)) and pcov is not None:
        s_sq = (criteria.fun(popt, *(ydata,stress))**2).sum()/(len(ydata)-len(initialguess))
        pcov = pcov * s_sq
      perr = np.sqrt(np.diag(pcov))
      fitResults.append(popt.tolist())
      fitErrors .append(perr.tolist())

      popt = np.concatenate((np.array(popt), np.repeat(options.exponent,nExponent)))
      perr = np.concatenate((np.array(perr), np.repeat(0.0,nExponent)))

      damask.util.croak('Needed {} function calls for fitting'.format(infodict['nfev']))
    except Exception as detail:
      damask.util.croak(detail)
      pass
    return popt


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
  global myLoad
  loadNo=loadcaseNo()
  if not os.path.isfile('%s.load'%loadNo):
    damask.util.croak('Generating load case for simulation %s (%s)'%(loadNo,thread))
    f=open('%s.load'%loadNo,'w')
    f.write(myLoad.getLoadcase(loadNo))
    f.close()
    s.release()
  else: s.release()

  s.acquire()
  if not os.path.isfile('%s_%i.spectralOut'%(options.geometry,loadNo)):
    damask.util.croak('Starting simulation %i (%s)'%(loadNo,thread))
    s.release()
    damask.util.execute('DAMASK_spectral -g %s -l %i'%(options.geometry,loadNo))
  else: s.release()
   
  s.acquire()
  if not os.path.isfile('./postProc/%s_%i.txt'%(options.geometry,loadNo)):
    damask.util.croak('Starting post processing for simulation %i (%s)'%(loadNo,thread))
    s.release()
    try:
      damask.util.execute('postResults --cr f,p --co totalshear %s_%i.spectralOut'%(options.geometry,loadNo))
    except:
      damask.util.execute('postResults --cr f,p %s_%i.spectralOut'%(options.geometry,loadNo))
    damask.util.execute('addCauchy ./postProc/%s_%i.txt'%(options.geometry,loadNo))
    damask.util.execute('addStrainTensors -0 -v ./postProc/%s_%i.txt'%(options.geometry,loadNo))
    damask.util.execute('addMises -s Cauchy -e ln(V) ./postProc/%s_%i.txt'%(options.geometry,loadNo))
  else: s.release()

  s.acquire()
  damask.util.croak('Reading values from simulation %i (%s)'%(loadNo,thread))
  s.release()

  refFile = './postProc/%s_%i.txt'%(options.geometry,loadNo)
  table = damask.ASCIItable(refFile)
  table.head_read()
  if options.fitting =='equivalentStrain':
    thresholdKey = 'Mises(ln(V))'
  elif options.fitting =='totalshear':
    thresholdKey = 'totalshear'
  s.acquire()
  for l in [thresholdKey,'1_Cauchy']:
    if l not in table.labels: damask.util.croak('%s not found'%l)
  s.release()
  table.data_readArray(['%i_Cauchy'%(i+1) for i in xrange(9)]+[thresholdKey]+['%i_ln(V)'%(i+1) for i in xrange(9)])

  line = 0
  lines = np.shape(table.data)[0]
  validity = np.zeros((int(options.yieldValue[2])), dtype=bool)
  yieldStress     = np.empty((int(options.yieldValue[2]),6),'d')
  deformationRate = np.empty((int(options.yieldValue[2]),6),'d')
  for i,threshold in enumerate(np.linspace(options.yieldValue[0],options.yieldValue[1],options.yieldValue[2])):
    while line < lines:
      if abs(table.data[line,9])>= threshold:
        upper,lower = abs(table.data[line,9]),abs(table.data[line-1,9])                             # values for linear interpolation
        stress = np.array(table.data[line-1,0:9] * (upper-threshold)/(upper-lower) + \
                          table.data[line  ,0:9] * (threshold-lower)/(upper-lower)).reshape(3,3)    # linear interpolation of stress values
        #stress = 0.5*(stress+stress.T)                                                             # symmetrise
        #for the mapping, a fuction from DAMASK (33to6) simplifies
        dstrain= np.array(table.data[line,10:] - table.data[line-1,10:]).reshape(3,3)
#
        yieldStress[i,0]= stress[0,0]; yieldStress[i,1]=stress[1,1]; yieldStress[i,2]=stress[2,2]
        yieldStress[i,3]=(stress[0,1] + stress[1,0])/2.0     #   0  3  5
        yieldStress[i,4]=(stress[1,2] + stress[2,1])/2.0     #   *  1  4  yieldStress
        yieldStress[i,5]=(stress[2,0] + stress[0,2])/2.0     #   *  *  2

#       D*dt = 0.5(L+L^T)*dt = 0.5*d(lnF + lnF^T) = dlnV
        deformationRate[i,0]= dstrain[0,0]; deformationRate[i,1]=dstrain[1,1]; deformationRate[i,2]=dstrain[2,2]
        deformationRate[i,3]=(dstrain[0,1] + dstrain[1,0])/2.0     #   0  3  5
        deformationRate[i,4]=(dstrain[1,2] + dstrain[2,1])/2.0     #   *  1  4
        deformationRate[i,5]=(dstrain[2,0] + dstrain[0,2])/2.0     #   *  *  2
        validity[i] = True
        break
      else:
        line+=1
    if not validity[i]:
      damask.util.croak('The data of result %i at the threshold %f is invalid,'\
                        +'the fitting at this point is skipped'%(loadNo,threshold))
  s.acquire()
  global stressAll, strainAll
  f=open(options.geometry+'_'+options.criterion+'.txt','w')
  f.write(' '.join([options.fitting]+myFit.report_labels())+'\n')
  try:
    for i,threshold in enumerate(np.linspace(options.yieldValue[0],options.yieldValue[1],options.yieldValue[2])):
      if validity[i]:
        a = (yieldStress[0][2]/stressUnit)**2 + (yieldStress[0][4]/stressUnit)**2 + (yieldStress[0][5]/stressUnit)**2
        stressAll[i]=np.append(stressAll[i], yieldStress[i]/stressUnit)
        strainAll[i]=np.append(strainAll[i], deformationRate[i])
        f.write( str(threshold)+' '+
                 ' '.join(map(str,myFit.fit(stressAll[i].reshape(len(stressAll[i])//6,6).transpose())))+'\n')
  except Exception as detail:
    damask.util.croak('Could not fit results of simulation (%s)'%thread)
    s.release()
    return
  damask.util.croak('\n')
  s.release()

def loadcaseNo():
  global N_simulations
  N_simulations+=1
  return N_simulations

def converged():
  global N_simulations; fitResidual

  if N_simulations < options.max:
    if len(fitResidual) > 5:
      residualList = np.array(fitResidual[len(fitResidual)-5:])
      if np.std(residualList)/np.max(residualList) < 0.05: 
        return True
    return False
  else:
    return True

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Performs calculations with various loads on given geometry file and fits yield surface.

""", version = scriptID)

# maybe make an option to specifiy if 2D/3D fitting should be done?

parser.add_option('-l','--load' ,      dest='load', type='float', nargs=3,
                                       help='load: final strain; increments; time %default', metavar='float int float')
parser.add_option('-g','--geometry',   dest='geometry', type='string',
                                       help='name of the geometry file [%default]', metavar='string')
parser.add_option('-c','--criterion',  dest='criterion', choices=fitCriteria.keys(),
                                       help='criterion for stopping simulations [%default]', metavar='string')
# best/worse fitting? Stopping?
parser.add_option('-f','--fitting',    dest='fitting', choices=thresholdParameter,
                                       help='yield criterion [%default]', metavar='string')
parser.add_option('-y','--yieldvalue', dest='yieldValue', type='float', nargs=3,
                                       help='yield points: start; end; count %default', metavar='float float int')
parser.add_option('--min',             dest='min', type='int',
                                       help='minimum number of simulations [%default]', metavar='int')
parser.add_option('--max',             dest='max', type='int',
                                       help='maximum number of iterations [%default]',  metavar='int')
parser.add_option('-t','--threads',    dest='threads', type='int',
                                       help='number of parallel executions [%default]',  metavar='int')
parser.add_option('-b','--bound',      dest='bounds', type='float', nargs=2,
                                       help='yield points: start; end; count %default', metavar='float float')
parser.add_option('-d','--dimension',  dest='dimension', type='int',
                                       help='dimension of the virtual test [%default]',  metavar='int')
parser.add_option('-v', '--vegter',    dest='vegter', action='store_true',
                                       help='Vegter criteria [%default]',  metavar='float') 
parser.add_option('-e', '--exponent',  dest='exponent', type='float',
                                       help='exponent of non-quadratic criteria',  metavar='int') 
parser.add_option('-u', '--uniaxial',  dest='eqStress', type='float',
                                       help='Equivalent stress',  metavar='float') 

parser.set_defaults(min        = 12)
parser.set_defaults(max        = 30)
parser.set_defaults(threads    = 4)
parser.set_defaults(yieldValue = (0.002,0.004,2))
parser.set_defaults(load       = (0.010,100,100.0))
parser.set_defaults(criterion  = 'vonmises')
parser.set_defaults(fitting    = 'totalshear')
parser.set_defaults(geometry   = '20grains16x16x16')
parser.set_defaults(bounds     = None)
parser.set_defaults(dimension  = 3)
parser.set_defaults(vegter     = 'False')
parser.set_defaults(exponent   = -1.0)

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
if options.dimension not in [2,3]:
  parser.error('Dimension is wrong, should be 2(plane stress state) or 3(general stress state)')
if not os.path.isfile('numerics.config'):
  damask.util.croak('numerics.config file not found')

if not os.path.isfile('material.config'):
  damask.util.croak('material.config file not found')


if options.vegter is True: 
  options.dimension = 2

if options.criterion == 'hill1948': stressUnit = 1.0e9
else                              : stressUnit = 1.0e6
fit_allResults = {}
if options.criterion == 'all':
  noExpoList = ['tresca', 'vonmises', 'hill1948', 'drucker']
  for criter in noExpoList:
    run = runFit(options.exponent, 0.0, fitCriteria[criter]['dimen'], criter)
    fit_allResults[criter] = {'results': fitResults, 'errors': fitErrors, 'residual': fitResidual}
  for criter in [x for x in fitCriteria.keys() if x not in noExpoList+['vegter','all'] ]:
    if options.eqStress == None:  
      # the user do not provide the equivalent stress, the equivalent stress will be fitted by von Mises
      eqStress = fit_allResults['vonmises']['results'][-1]
    run = runFit(options.exponent, eqStress, fitCriteria[criter]['dimen'], criter)
    fit_allResults[criter] = {'results': fitResults, 'errors': fitErrors, 'residual': fitResidual}
else:
  criter = options.criterion
  if fitCriteria[criter]['nExpo'] == 0:
    run = runFit(options.exponent, 0.0, fitCriteria[criter]['dimen'], criter)
  else:
    if options.eqStress == None: 
      # the user do not provide the equivalent stress, the equivalent stress will be fitted by von Mises
      run = runFit(options.exponent, 0.0, fitCriteria[criter]['dimen'], 'vonmises')
      eqStress = fitResults[-1]
    run = runFit(options.exponent, eqStress, fitCriteria[criter]['dimen'], criter)
  fit_allResults[criter] = {'results': fitResults, 'errors': fitErrors, 'residual': fitResidual}

damask.util.croak('Finished fitting to yield criteria')
