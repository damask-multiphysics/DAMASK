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

def principalStress(p):
  sin = np.sin;    cos = np.cos
  I1,I2,I3 = invariant(p)

  third = 1.0/3.0
  I1s3I2= (I1**2 - 3.0*I2)**0.5
  numer = 2.0*I1**3 - 9.0*I1*I2 + 27.0*I3
  denom = I1s3I2**(-3.0)
  cs    = 0.5*numer*denom
  phi   = np.arccos(cs)/3.0
  t1    = I1/3.0; t2 = 2.0/3.0*I1s3I2
  return np.array( [t1 + t2*cos(phi), t1+t2*cos(phi+np.pi*2.0/3.0), t1+t2*cos(phi+np.pi*4.0/3.0)])

def principalStrs_Der(p, (s1, s2, s3, s4, s5, s6), Karafillis=False):
  sin = np.sin;    cos = np.cos
  I1,I2,I3 = invariant(p)

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
  one    = np.ones_like(s1); zero = np.zeros_like(s1); dim  = len(s1)
  dIdp = np.array([[one,    one,    one,    zero,   zero,   zero], 
                   [p[1]+p[2], p[2]+p[0], p[0]+p[1], -2.0*p[3], -2.0*p[4], -2.0*p[5]],
                   [p[1]*p[2]-p[4]**2, p[2]*p[0]-p[5]**2, p[0]*p[1]-p[3]**2, 
                    -2.0*p[3]*p[2]+2.0*p[4]*p[5], -2.0*p[4]*p[0]+2.0*p[5]*p[3], -2.0*p[5]*p[1]+2.0*p[3]*p[4]] ]) 
  if Karafillis:
    dpdc = np.array([[zero,s2-s3,s3-s2], [s1-s3,zero,s3-s1], [s1-s2,s2-s1,zero]])
    dSdp = np.array([np.dot(dSdI[:,:,i],dIdp[:,:,i]).T for i in xrange(dim)]).T
    return np.concatenate((np.array([np.dot(dSdp[:,0:3,i], dpdc[:,:,i].T).T/3.0 for i in xrange(dim)]).T,
                           np.vstack([dSdp[:,3]*s4,dSdp[:,4]*s5,dSdp[:,5]*s6]).T.reshape(dim,3,3).T), axis=1)
  else:
    dIdc=np.array([[-dIdp[i,0]*s2, -dIdp[i,1]*s1, -dIdp[i,1]*s3,
                    -dIdp[i,2]*s2, -dIdp[i,2]*s1, -dIdp[i,0]*s3,
                     dIdp[i,3]*s4,  dIdp[i,4]*s5,  dIdp[i,5]*s6 ] for i in xrange(3)])
    return np.array([np.dot(dSdI[:,:,i],dIdc[:,:,i]).T for i in xrange(dim)]).T

def invariant(sigmas):
    s11,s22,s33,s12,s23,s31 = sigmas
    I1  = s11 + s22 + s33
    I2  = s11*s22 + s22*s33 + s33*s11 - s12**2 - s23**2 - s31**2
    I3  = s11*s22*s33 + 2.0*s12*s23*s31 - s12**2*s33 - s23**2*s11 - s31**2*s22
    return (I1,I2,I3)

def formatOutput(n, type='%-14.6f'):
  return ''.join([type for i in xrange(n)])

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

class Criteria(object):
  '''
    residuum of anisotropic Barlat 1991 yield criterion (eq. 2.37)
  '''
  def __init__(self, criterion, uniaxialStress,exponent):
    self.stress0 = uniaxialStress
    if exponent < 0.0:
      self.mFix = [False, exponent]
    else:
      self.mFix = [True, exponent]
    self.func     = fitCriteria[criterion]['func']
    self.criteria = criterion
  def fun(self, paras, ydata, sigmas):
    return self.func(self.stress0, paras, sigmas,self.mFix,self.criteria)
  def jac(self, paras, ydata, sigmas):
    return self.func(self.stress0, paras, sigmas,self.mFix,self.criteria,Jac=True)

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

def Tresca(eqStress, paras, sigmas, mFix, criteria, Jac = False):
  '''
    Tresca yield criterion
    the fitted parameters is: paras(sigma0)
  '''
  lambdas = principalStresses(sigmas)
  r = np.amax(np.array([abs(lambdas[2,:]-lambdas[1,:]),\
                        abs(lambdas[1,:]-lambdas[0,:]),\
                        abs(lambdas[0,:]-lambdas[2,:])]),0) - paras
  if not Jac:
    return r.ravel()
  else:
    return -np.ones(len(r))

def Cazacu_Barlat(eqStress, paras, sigmas, mFix, criteria, Jac = False):
  '''
    Cazacu–Barlat (CB) yield criterion
    the fitted parameters are:
      a1,a2,a3,a6; b1,b2,b3,b4,b5,b10; c for plane stress
      a1,a2,a3,a4,a5,a6; b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11; c: for general case
    mFix are invalid input
  '''
  if criteria == 'cb2d':
    coeffa, coeffb, c = paras[0:4],paras[4:10],paras[10]
  else:
    coeffa, coeffb, c = paras[0:6],paras[6:17],paras[17]

  s11,s22,s33,s12,s23,s31 = sigmas
  if criteria == 'cb2d': s33=s23=s31 = np.zeros_like(s11)
  s1_2, s2_2, s3_2, s12_2, s23_2, s31_2 = np.array([s11,s22,s33,s12,s23,s31])**2
  s1_3, s2_3, s3_3, s123, s321 = s11*s1_2, s22*s2_2, s33*s3_2,s11*s22*s33, s12*s23*s31
  d12,d23,d31 = s11-s22, s22-s33, s33-s11
  
  jb1 = (s1_3 + 2.0*s3_3)/27.0 - s22*s1_2/9.0 - (s11+s22)*s3_2/9.0 + s123/4.5
  jb2 = (s1_3 -     s3_3)/27.0 - s33*s1_2/9.0 +  s11     *s3_2/9.0
  jb3 = (s2_3 -     s3_3)/27.0 - s33*s2_2/9.0 +  s22     *s3_2/9.0
  jb4 = (s2_3 + 2.0*s3_3)/27.0 - s11*s2_2/9.0 - (s11+s22)*s3_2/9.0 + s123/4.5

  jb5, jb10 = -d12*s12_2/3.0, -d31*s12_2/1.5
  jb6, jb7  = -d12*s23_2/3.0,  d31*s23_2/3.0
  jb8, jb9  =  d31*s31_2/3.0,  d12*s31_2/1.5
  jb11      =  s321*2.0
  
  if criteria == 'cb3d':
    dJ2da = np.array([d12**2/6.0, d23**2/6.0, d31**2/6.0, s12_2,s23_2,s31_2])
    dJ3db = np.array([jb1,jb2,jb3,jb4,jb5,jb6,jb7,jb8,jb9,jb10,jb11])
  else:  # plane stress
    dJ2da = np.array([d12**2/6.0, s2_2/6.0, s1_2/6.0, s12_2])
    dJ3db = np.array([jb1,jb2,jb3,jb4,jb5,jb10])

  J20 = np.dot(coeffa,dJ2da)
  J30 = np.dot(coeffb,dJ3db)
  f0  = (J20**3 - c*J30**2)/18.0
  r   = f0**(1.0/6.0)*(3.0/eqStress)

  if not Jac:
    return (r - 1.0).ravel()
  else:
    df = r/f0/108.0
    return np.vstack((df*3.0*J20**2.0*dJ2da, -df*2.0*J30*c*dJ3db, -df*J30**2)).T

def Drucker(eqStress, paras, sigmas, mFix, criteria, Jac = False):
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
  I1,I2,I3 = invariant(sigmas)
  J2  = I1**2/3.0 - I2
  J3  = I1**3/13.5 - I1*I2/3.0 + I3
  J2_3p = J2**(3.0*p)      
  J3_2p = J3**(2.0*p)
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

def Hill1948(eqStress, paras, sigmas, mFix, criteria, Jac = False):
  '''
    Hill 1948 yield criterion
    the fitted parameters are F, G, H, L, M, N
    eqStress, criteria are invalid input
  '''
  s11,s22,s33,s12,s23,s31 = sigmas
  jac = np.array([(s22-s33)**2,(s33-s11)**2,(s11-s22)**2, 2.0*s23**2,2.0*s31**2,2.0*s12**2])
  if not Jac:
    return (np.dot(paras,jac)/2.0-0.5).ravel()
  else:
    return jac.T

def Hill1979(eqStress,paras, sigmas, mFix, criteria, Jac = False):
  '''
    Hill 1979 yield criterion
    the fitted parameters are: f,g,h,a,b,c,m
  '''
  if mFix[0]: m = mFix[1]
  else:       m = paras[-1]

  coeff  = paras[0:6]
  s1,s2,s3 = principalStresses(sigmas)
  diffs  = np.array([s2-s3, s3-s1, s1-s2, 2.0*s1-s2-s3, 2.0*s2-s3-s1, 2.0*s3-s1-s2])**2 
  diffsm = diffs**(m/2.0)
  base   = np.dot(coeff,diffsm)
  r = base**(1.0/m)/eqStress #left = base**mi

  if not Jac:
    return (r-1.0).ravel()
  else:
    drdb = r/base/m
    dbdm = np.dot(coeff,diffsm*math_ln(diffs)) #****0.5
    jm = drdb*dbdm + r*math_ln(base)/(-m**2)

    if mFix[0]: return np.vstack((drdb*diffsm)).T
    else:       return np.vstack((drdb*diffsm, jm)).T

def Hosford(eqStress, paras, sigmas, mFix, criteria, Jac = False):
  '''
    Hosford family criteria
    the fitted parameters are:
      von Mises: sigma0
      Hershey: (1) sigma0, a, when a is not fixed; (2) sigma0, when a is fixed
      general Hosford: (1) F,G,H, a, when a is not fixed; (2) F,G,H, when a is fixed
  '''

  if criteria == 'vonmises':
    coeff = np.ones(3)
    a = 2.0
    sigma0 = paras
  elif criteria == 'hershey':
    sigma0 = paras[0]
    if mFix[0]: a = mFix[1]
    else:       a = paras[1]
    coeff = np.ones(3)
  else:
    print '11'
    coeff = paras[0:3]
    if mFix[0]: a = mFix[1]
    else:       a = paras[3]

  s1,s2,s3 = principalStresses(sigmas)
  diffs  = np.abs(np.array([s2-s3, s3-s1, s1-s2]))
  diffsm = diffs**a
  base   = np.dot(coeff,diffsm)
  expo   = 1.0/a
  r      = (base/2.0)**expo/sigma0

  if not Jac:
    return (r-1.0).ravel()
  else:
    if criteria == 'vonmises': # von Mises
      return -r/sigma0
    else:
      dbda = np.dot(coeff,diffsm*math_ln(diffs))
      drdb = r/base*expo
      ja = drdb*dbda + r*math_ln(base/2.0)*(-expo*expo)
      if criteria == 'hershey':  # Hershey
        if mFix[0]: return -r/sigma0
        else:       return np.vstack((-r/sigma0, ja)).T
      else:                                        # Anisotropic Hosford
        if mFix[0]: return np.vstack((drdb*diffsm)).T
        else:       return np.vstack((drdb*diffsm, ja)).T

def Barlat1991(eqStress, paras, sigmas, mFix, criteria, Jac=False):
  '''
    Barlat 1991 criteria
    the fitted parameters are:
      Isotropic: sigma0, m
      Anisotropic: a, b, c, f, g, h, m
  '''
  if criteria == 'barlat1991iso':
    sigma0 = paras[0]
    coeff  = np.ones(6)
  else:
    sigma0 = eqStress
    coeff  = paras[0:6]   
  if mFix[0]:m = mFix[1]
  else:      m = paras[-1]
      
  cos = np.cos; sin = np.sin; pi = np.pi; abs = np.abs
  s1,s2,s3,s4,s5,s6 = sigmas
  dXdx = np.array([s2-s3,s3-s1,s1-s2,s5,s6,s4])
  A,B,C,F,G,H = np.array(coeff)[:,None]*dXdx

  I2 = (F*F + G*G + H*H)/3.0+ ((A-C)**2+(C-B)**2+(B-A)**2)/54.0
  I3 = (C-B)*(A-C)*(B-A)/54.0 + F*G*H - ((C-B)*F*F + (A-C)*G*G + (B-A)*H*H)/6.0
  phi1 = np.arccos(I3/I2**1.5)/3.0 + pi/6.0; absc1 = 2.0*abs(cos(phi1))
  phi2 = phi1                      + pi/3.0; absc2 = 2.0*abs(cos(phi2))
  phi3 = phi2                      + pi/3.0; absc3 = 2.0*abs(cos(phi3))
  left = ( absc1**m + absc2**m + absc3**m )/2.0
  r    = left**(1.0/m)*np.sqrt(3.0*I2)/sigma0

  if not Jac:
    return (r - 1.0).ravel()
  else:
    dfdl = r/left/m
    jm = r*math_ln(left)/(-m**2) + dfdl*0.5*(
         absc1**m*math_ln(absc1) + absc2**m*math_ln(absc2) + absc3**m*math_ln(absc3) ) 
    if criteria == 'barlat1991iso':
      js = -(r + 1.0)/sigma0
      if mFix[0]: return js
      else:       return np.vstack((js,jm)).T
    else:
      da,db,dc = (2.0*A-B-C)/18.0, (2.0*B-C-A)/18.0, (2.0*C-A-B)/18.0
      dI2dx = np.array([da, db, dc, F,G,H])/1.5*dXdx
      dI3dx = np.array([da*(B-C) + (H**2-G**2)/2.0, db*(C-A) + (F**2-H**2)/2.0, dc*(A-B) + (G**2-F**2)/2.0,
                        (H*G + (B-C))*F, (F*H + (C-A))*G, (G*F + (A-B))*H])/3.0*dXdx
      darccos = -(1.0 - I3**2/I2**3)**(-0.5)

      dfdc  = dfdl*0.5*m
      dfdcos = lambda phi : dfdc*(2.0*abs(cos(phi)))**(1.0/m-1.0)*np.sign(cos(phi))*(-sin(phi)/1.5)
      dfdthe= (dfdcos(phi1) + dfdcos(phi2) + dfdcos(phi3)) 
      dfdI2 = dfdthe*darccos*I3*(-1.5)*I2**(-2.5);   dfdI3 = dfdthe*darccos*I2**(-1.5)

      if mFix[0]: return np.vstack((dfdI2*dI2dx+dfdI3*dI3dx)).T
      else:       return np.vstack((dfdI2*dI2dx+dfdI3*dI3dx, jm)).T

def BBC2003(eqStress, paras, sigmas, mFix, criteria, Jac=False):
  '''
    residuum of the BBC2003 yield criterion for plain stress
  '''
  a,b,c, d,e,f,g= paras[0:7]
  if mFix[0]: k = mFix[1]
  else:       k = paras[-1]

  s11 = sigmas[0]; s22 = sigmas[1]; s12 = sigmas[3]
  k2  = 2.0*k
  M = d+e;  N = e+f;  P = (d-e)/2.0;  Q = (e-f)/2.0; R = g**2
  Gamma =    M*s11 + N*s22
  Psi   = ( (P*s11 + Q*s22)**2 + s12**2*R )**0.5

  l1  = b*Gamma + c*Psi;   l1s = l1**2
  l2  = b*Gamma - c*Psi;   l2s = l2**2
  l3  = 2.0*c*Psi;         l3s = l3**2
                       
  left = a*l1s**k + a*l2s**k + (1-a)*l3s**k
  sBar = left**(1.0/k2);  r = sBar/eqStress - 1.0
  if not Jac:
    return r.ravel()
  else:
    temp = (P*s11 + Q*s22)/Psi
    dPsidP = temp*s11;  dPsidQ = temp*s22;  dPsidR = 0.5*s12**2/Psi
    expo = 0.5/k;  k1 = k-1.0

    dsBardl = expo*sBar/left/eqStress
    dsBarde = sBar*math_ln(left);   dedk = expo/(-k)
    dldl1 =    a *k*(l1s**k1)*(2.0*l1)
    dldl2 =    a *k*(l2s**k1)*(2.0*l2)
    dldl3 = (1-a)*k*(l3s**k1)*(2.0*l3)

    dldGama = (dldl1 + dldl2)*b
    dldPsi  = (dldl1 - dldl2 + 2.0*dldl3)*c

    dlda = l1s**k + l2s**k - l3s**k
    dldb = dldl1*Gamma + dldl2*Gamma
    dldc = dldl1*Psi   - dldl2*Psi + dldl3*2.0*Psi
    dldk = a*math_ln(l1s)*l1s**k + a*math_ln(l2s)*l2s**k + (1-a)*math_ln(l3s)*l3s**k

    ja = dsBardl * dlda
    jb = dsBardl * dldb
    jc = dsBardl * dldc
    jd = dsBardl *(dldGama*s11 + dldPsi*dPsidP*0.5)
    je = dsBardl *(dldGama*(s11+s22) + dldPsi*(dPsidP*(-0.5) + dPsidQ*0.5) )
    jf = dsBardl *(dldGama*s22 + dldPsi*dPsidQ*(-0.5))
    jg = dsBardl * dldPsi * dPsidR * 2.0*g
    jk = dsBardl * dldk + dsBarde * dedk

    if mFix[0]: return np.vstack((ja,jb,jc,jd, je, jf,jg)).T
    else:       return np.vstack((ja,jb,jc,jd, je, jf,jg,jk)).T

def BBC2005(eqStress, paras, sigmas, mFix, criteria, Jac=False):
  '''
    residuum of the BBC2005 yield criterion for plain stress
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

def Yld200418p(eqStress, paras, sigmas, mFix, criteria, Jac=False):
  '''
    C: c12,c21,c23,c32,c13,c31,c44,c55,c66
    D: d12,d21,d23,d32,d31,d13,d44,d55,d66
  '''
  C,D = paras[0:9], paras[9:18]
  if mFix[0]: m = mFix[1]
  else:       m = paras[-1]

  sv = (sigmas[0] + sigmas[1] + sigmas[2])/3.0
  sdev = np.vstack((sigmas[0:3]-sv,sigmas[3:6]))
  ys = lambda sdev, C: np.array([-C[0]*sdev[1]-C[5]*sdev[2], -C[1]*sdev[0]-C[2]*sdev[2], 
                                 -C[4]*sdev[0]-C[3]*sdev[1],  C[6]*sdev[3],C[7]*sdev[4], C[8]*sdev[5]])
  p,q = ys(sdev, C), ys(sdev, D)
  pLambdas, qLambdas = principalStress(p), principalStress(q)   # no sort

  m2 = m/2.0;  m1 = 1.0/m;  m21 = m2-1.0; x3 = xrange(3); dim = len(sv)
  PiQj  = np.array([(pLambdas[i,:]-qLambdas[j,:]) for i in x3 for j in x3])
  QiPj  = np.array([(qLambdas[i,:]-pLambdas[j,:]) for i in x3 for j in x3]).reshape(3,3,dim)
  PiQjs = PiQj**2
  phi   = np.sum(PiQjs**m2,axis=0)
  r  = (0.25*phi)**m1/eqStress

  if not Jac:
    return (r - 1.0).ravel()
  else:
    drdphi =  r*m1/phi*4.0
    dphidm = np.sum(PiQjs**m2*math_ln(PiQjs),axis=0)*0.5
    dPdc, dQdd = principalStrs_Der(p, sdev), principalStrs_Der(q, sdev)
    PiQjs3d = (PiQjs**m21).reshape(3,3,dim)
    dphidP = -m*np.array([np.diag(np.dot(PiQjs3d[:,:,i], QiPj   [:,:,i])) for i in xrange(dim)]).T
    dphidQ =  m*np.array([np.diag(np.dot(QiPj   [:,:,i], PiQjs3d[:,:,i])) for i in xrange(dim)]).T

    jm = drdphi*dphidm + r*math_ln(0.25*phi)*(-m1*m1)
    jc = drdphi*np.sum([dphidP[i]*dPdc[i] for i in x3],axis=0)
    jd = drdphi*np.sum([dphidQ[i]*dQdd[i] for i in x3],axis=0)

    if mFix[0]: return np.vstack((jc,jd)).T
    else:       return np.vstack((jc,jd,jm)).T

def KarafillisBoyce(eqStress, paras, sigmas, mFix, criteria, Jac=False):
  ks = lambda (s1,s2,s3,s4,s5,s6),(c1,c2,c3,c4,c5,c6): np.array( [
         ((c2+c3)*s1-c3*s2-c2*s3)/3.0,  ((c3+c1)*s2-c3*s1-c1*s3)/3.0,
         ((c1+c2)*s3-c2*s1-c1*s2)/3.0,  c4*s4, c5*s5, c6*s6 ])

  C1,C2,alpha  = paras[0:6], paras[6:12], paras[12]
  if mFix[0]: b1=b2=a = mFix[1]
  else:       b1,b2,a = paras[12:15]
  print b1,b2,a
  p,q = ks(sigmas, C1), ks(sigmas, C2)
  plambdas,qlambdas = principalStress(p), principalStress(q)
  b1i,b2i,ai,rb2 = 1.0/b1, 1.0/b2, 1.0/a, 3.0**b2/(2.0**b2+2.0)

  difP  = np.array([plambdas[1]-plambdas[2], plambdas[2]-plambdas[0], plambdas[0]-plambdas[1]])
  difPs = difP**2; difPb1 = difPs**(b1/2.0-1.0)
  Qs    = qlambdas**2

  phi10, phi20 = np.sum(difPs**(b1/2.0),axis = 0), np.sum(Qs**(b2/2.0),axis = 0)
  phi1,  phi2  = (0.5*phi10)**b1i, (rb2*phi20)**b2i
  Stress = alpha*phi1**a + (1.0-alpha)*phi2**a
  r      = Stress**ai/eqStress

  if not Jac:
    return (r-1.0).ravel()
  else:
    drds = r*ai/Stress
    dsda = alpha*phi1**a*math_ln(phi1) + (1.0-alpha)*phi2**a*math_ln(phi2)

    dphi1dP = phi1/phi10*np.array([       -difPb1[1]*difP[1]+difPb1[2]*difP[2],
      difPb1[0]*difP[0]-difPb1[2]*difP[2], difPb1[1]*difP[1]-difPb1[0]*difP[0]])
    dphi2dQ = phi2/phi20*Qs*qlambdas*(b2/2.0-1.0)
    dPdc = principalStrs_Der(p, sigmas, Karafillis=True)
    dQdc = principalStrs_Der(q, sigmas, Karafillis=True)
    dphi10db1 = np.sum(difPs**(b1/2.0)*math_ln(difPs), axis=0)*0.5
    dphi20db2 = np.sum(   Qs**(b2/2.0)*math_ln(   Qs), axis=0)*0.5

    drb2db2   = rb2*math_ln(3.0) - rb2*math_ln(2.0)/(1.0+2.0**(1.0-b2))
    dphi1db1  = phi1*math_ln(phi10)*(-b1i*b1i) + b1i*phi1/(0.5*phi10)* 0.5*dphi10db1
    dphi2db2  = phi2*math_ln(phi20)*(-b2i*b2i) + b2i*phi2/(rb2*phi20)*(rb2*dphi20db2 + drb2db2*phi20)
    ja  = drds*dsda - r*math_ln(Stress)/a/a  #drda
    jb1 = dphi1db1*(drds*a*phi1**(a-1)*alpha ) 
    jb2 = dphi2db2*(drds*a*phi2**(a-1)*(1.0-alpha))
    jc1 = np.sum([dphi1dP[i]*dPdc[i] for i in xrange(3)],axis=0)*drds*a*phi1**(a-1.0)*alpha 
    jc2 = np.sum([dphi2dQ[i]*dQdc[i] for i in xrange(3)],axis=0)*drds*a*phi2**(a-1.0)*(1.0-alpha)
    jalpha = drds * (phi1**a - phi2**a)

    if mFix[0]: return np.vstack((jc1,jc2,jalpha)).T
    else:       return np.vstack((jc1,jc2,jalpha,jb1,jb2,ja)).T


fitCriteria = {
  'tresca'         :{'func' : Tresca,
                     'nExpo': 0,'err':np.inf,
                     'bound': [(None,None)],
                     'paras': 'Initial yield stress:',
                     'text' : '\nCoefficient of Tresca criterion:\nsigma0: ',
                     'error': 'The standard deviation error is: '
                    },
  'vonmises'       :{'func' : Hosford,
                     'nExpo': 0,'err':np.inf,
                     'bound': [(None,None)],
                     'paras': 'Initial yield stress:',
                     'text' : '\nCoefficient of Huber-Mises-Hencky criterion:\nsigma0: ',
                     'error': 'The standard deviation error is: '
                    },
  'hershey'        :{'func' : Hosford,
                     'nExpo': 1,'err':np.inf,
                     'bound': [(None,None)]+[(2.0,8.0)],
                     'paras': 'Initial yield stress, a:',
                     'text' : '\nCoefficients of Hershey criterion:\nsigma0, a: ',
                     'error': 'The standard deviation errors are: '
                    },
  'ghosford'       :{'func' : Hosford,
                     'nExpo': 1,'err':np.inf,
                     'bound': [(0.0,2.0)]*3+[(0.0,12.0)],
                     'paras': 'F, G, H, a:',
                     'text' : '\nCoefficients of Hosford criterion:F, G, H, a: ',
                     'error': 'The standard deviation errors are: '
                    },
  'hill1948'       :{'func' : Hill1948,
                     'nExpo': 0,'err':np.inf,
                     'bound': [(None,None)]*6,
                     'paras': 'Normalized [F, G, H, L, M, N]:',
                     'text' : '\nCoefficients of Hill1948 criterion:\n[F, G, H, L, M, N]:'+' '*16,
                     'error': 'The standard deviation errors are: '
                    },
  'hill1979'       :{'func' : Hill1979,
                     'nExpo': 1,'err':np.inf,
                     'bound': [(-2.0,2.0)]*6+[(1.0,8.0)],
                     'paras': 'f,g,h,a,b,c,m:',
                     'text' : '\nCoefficients of Hill1979 criterion:\n f,g,h,a,b,c,m:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'drucker'        :{'func' : Drucker,
                     'nExpo': 0,'err':np.inf,
                     'bound': [(None,None)]*2,
                     'paras': 'Initial yield stress, C_D:',
                     'text' : '\nCoefficients of Drucker criterion:\nsigma0, C_D: ',
                     'error': 'The standard deviation errors are: '
                    },
  'gdrucker'       :{'func' : Drucker,
                     'nExpo': 1,'err':np.inf,
                     'bound': [(None,None)]*2+[(0.0,6.0)],
                     'paras': 'Initial yield stress, C_D, p:',
                     'text' : '\nCoefficients of general Drucker criterion:\nsigma0, C_D, p: ',
                     'error': 'The standard deviation errors are: '
                    },
  'barlat1991iso'  :{'func' : Barlat1991,
                     'nExpo': 1,'err':np.inf,
                     'bound': [(None,None)]+[(0.0,12.0)],
                     'paras': 'Initial yield stress, m:',
                     'text' : '\nCoefficients of isotropic Barlat 1991 criterion:\nsigma0, m:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'barlat1991aniso':{'func' : Barlat1991,
                     'nExpo': 1,'err':np.inf,
                     'name' : 'Barlat1991',
                     'bound': [(None,None)]*6+[(0.0,12.0)],
                     'text' : '\nCoefficients of anisotropic Barlat 1991 criterion: a, b, c, f, g, h, m:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'bbc2003'        :{'func' : BBC2003,
                     'nExpo': 1,'err':np.inf,
                     'bound': [(None,None)]*7+[(1.0,8.0)],
                     'paras': 'a, b, c, d, e, f, g, k:',
                     'text' : '\nCoefficients of Banabic-Balan-Comsa 2003 criterion: a, b, c, d, e, f, g, k:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'bbc2005'        :{'func' : BBC2005,
                     'nExpo': 1,'err':np.inf,
                     'bound': [(None,None)]*8+[(0.0,12.0)],
                     'paras': 'a, b, L ,M, N, P, Q, R, k:',
                     'text' : '\nCoefficients of Banabic-Balan-Comsa 2005 criterion: a, b, L ,M, N, P, Q, R, k:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'cb2d'           :{'func' : Cazacu_Barlat,
                     'nExpo': 0,'err':np.inf,
                     'bound': [(None,None)]*11,
                     'paras': 'a1,a2,a3,a6; b1,b2,b3,b4,b5,b10; c:',
                     'text' : '\nCoefficients of Cazacu Barlat yield criterion for plane stress: \
                               \n a1,a2,a3,a6; b1,b2,b3,b4,b5,b10; c:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'cb3d'           :{'func' : Cazacu_Barlat,
                     'nExpo': 0,'err':np.inf,
                     'bound': [(None,None)]*18,
                     'paras': 'a1,a2,a3,a4,a5,a6; b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11; c:',
                     'text' : '\nCoefficients of Cazacu Barlat yield criterion: \
                               \n a1,a2,a3,a4,a5,a6; b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11; c\n',
                     'error': 'The standard deviation errors are: '
                    },
  'yld200418p'     :{'func' : Yld200418p,
                     'nExpo': 1,'err':np.inf,
                     'bound': [(None,None)]*18+[(0.0,12.0)],
                     'paras': 'c12,c21,c23,c32,c31,c13,c44,c55,c66,d12,d21,d23,d32,d31,d13,d44,d55,d66,m:',
                     'text' : '\nCoefficients of Yld2004-18p yield criterion: \
                               \n c12,c21,c23,c32,c31,c13,c44,c55,c66,d12,d21,d23,d32,d31,d13,d44,d55,d66,m\n',
                     'error': 'The standard deviation errors are: '
                    },
  'karafillis'     :{'func' : KarafillisBoyce,
                     'nExpo': 3,'err':np.inf,
                     'bound': [(None,None)]*12+[(0.0,1.0)]+[(0.0,12.0)]*3,
                     'paras': 'c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,alpha,b1,b2,a',
                     'text' : '\nCoefficients of Karafillis-Boyce yield criterion: \
                               \n c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,alpha,b1,b2,a\n',
                     'error': 'The standard deviation errors are: '
                    }
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
    global fitResults; fitErrors
    if options.exponent > 0.0: nExponent = nExpo
    else:                      nExponent = 0
    nameCriterion = self.name.lower()
    criteria      = Criteria(nameCriterion,self.uniaxial,self.expo)
    textParas     = fitCriteria[nameCriterion]['text'] + formatOutput(numParas+nExponent)
    textError     = fitCriteria[nameCriterion]['error']+ formatOutput(numParas+nExponent,'%-14.8f')+'\n'
    bounds        = fitCriteria[nameCriterion]['bound']  # Default bounds, no bound
    guess0        = Guess                                # Default initial guess, depends on bounds

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
      fitErrors .append(perr.tolist())

      popt = np.concatenate((np.array(popt), np.repeat(options.exponent,nExponent)))
      perr = np.concatenate((np.array(perr), np.repeat(0.0,nExponent)))

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
parser.add_option('-c','--criterion',     dest='criterion', choices=fitCriteria.keys(),
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
parser.add_option('-e', '--exponent',dest='exponent', type='float',
                                     help='exponent of non-quadratic criteria') 
parser.add_option('-u', '--uniaxial',dest='eqStress', type='float',
                                     help='Equivalent stress') 
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
if not os.path.isfile('numerics.config'):
  print('numerics.config file not found')

if not os.path.isfile('material.config'):
  print('material.config file not found')

numParas = len(fitCriteria[options.criterion]['bound'])
nExpo = fitCriteria[options.criterion]['nExpo']
Guess = []
if options.exponent > 0.0: 
  numParas = numParas-nExpo # User defines the exponents
  fitCriteria[options.criterion]['bound'] = fitCriteria[options.criterion]['bound'][:numParas]
for i in xrange(numParas):
  temp = fitCriteria[options.criterion]['bound'][i]
  if fitCriteria[options.criterion]['bound'][i] == (None,None):Guess.append(1.0)
  else: 
    g = (temp[0]+temp[1])/2.0
    if g == 0: g = temp[1]*0.5
    Guess.append(g)
print Guess
print fitCriteria[options.criterion]['bound']
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
