#!/usr/bin/python
# -*- coding: UTF-8 no BOM -*-

import threading,time,os,subprocess,shlex,string
import numpy as np
from scipy.optimize import curve_fit
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
  def fun(self, sigma0, ydata, sigmas):
    return HosfordBasis(sigma0, 1.0,1.0,1.0, 2.0, sigmas)
  def jac(self, sigma0, ydata, sigmas):
    return HosfordBasis(sigma0, 1.0,1.0,1.0, 2.0, sigmas, Jac=True, nParas=1)

class Drucker(object):
  '''
    residuum of Drucker yield criterion (eq. 2.41, F = sigma0)
  '''
  def fun(sigma0, C_D, ydata, sigmas):
    return DruckerBasis(sigma0, C_D, 1.0, ydata, sigmas)
  def jac(sigma0, C_D, ydata, sigmas):
    pass

class generalDrucker(object):
  '''
    residuum of general Drucker yield criterion (eq. 2.42, F = sigma0)
  '''
  def fun(sigma0, C_D, ydata, sigmas):
    return DruckerBasis(sigma0, C_D, p, ydata, sigmas)
  def jac(sigma0, C_D, ydata, sigmas):
    pass

class Hosford(object):
  '''
    residuum of Hershey yield criterion (eq. 2.43, Y = sigma0)
  '''
  def fun(self, (sigma0, a), ydata, sigmas):
    return HosfordBasis(sigma0, 1.0,1.0,1.0, a, sigmas)
  def jac(self, (sigma0, a), ydata, sigmas):
    return HosfordBasis(sigma0, 1.0,1.0,1.0, a, sigmas, Jac=True, nParas=2)


#more to do
# KarafillisAndBoyce

# ---------------------------------------------------------------------------------------------
# isotropic yield surfaces
# ---------------------------------------------------------------------------------------------

class Hill1948(object):
  '''
    residuum of Hill 1948 quadratic yield criterion (eq. 2.48)
  '''
  def fun(self, (F,G,H,L,M,N), ydata, sigmas):
    r = F*(sigmas[1]-sigmas[2])**2.0 + G*(sigmas[2]-sigmas[0])**2.0 + H*(sigmas[0]-sigmas[1])**2.0\
      + 2.0*L*sigmas[4]**2.0 + 2.0*M*sigmas[5]**2.0 + 2.0*N*sigmas[3]**2.0 - 1.0
    return r.ravel()/2.0
  def jac(self, (F,G,H,L,M,N), ydata, sigmas):
    pass

#more to do
# Hill 1979

# Hill 1990,1993 need special stresses to fit

class generalHosford(object):
  '''
    residuum of Hershey yield criterion (eq. 2.104, sigmas = sigma0)
  '''
  def fun(self, (sigma0, F, G, H, a), ydata, sigmas, nParas=5):
    return HosfordBasis(sigma0, F, G, H, a, sigmas)
  def jac(self, (sigma0, F, G, H, a), ydata, sigmas):
    return HosfordBasis(sigma0, F,G,H, a, sigmas, Jac=True, nParas=5)

class Barlat1991iso(object):
  '''
    residuum of isotropic Barlat 1991 yield criterion (eq. 2.37)
  '''
  def fun(self, (sigma0, m), ydata, sigmas):
    return Barlat1991Basis(sigma0, 1.0,1.0,1.0,1.0,1.0,1.0, m, sigmas)
  def jac(self, (sigma0, m), ydata, sigmas):
    pass

class Barlat1991aniso(object):
  '''
    residuum of anisotropic Barlat 1991 yield criterion (eq. 2.37)
  '''
  def fun(self, (sigma0, a,b,c,f,g,h, m), ydata, sigmas):
    return Barlat1991Basis(sigma0, a,b,c,f,g,h, m, sigmas)
  def jac(self, (sigma0, a,b,c,f,g,h, m), ydata, sigmas):
    pass

def Cazacu_Barlat3D(sigma0,a1,a2,a3,a4,a5,a6, b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11, c,
                    ydata, sigmas):
  '''
    residuum of the Cazacu–Barlat (CZ) yield criterion
  '''
  s11 = sigmas[0]; s22 = sigmas[1]; s33 = sigmas[2]
  s12 = sigmas[3]; s23 = sigmas[4]; s31 = sigmas[5]
  
  J20 = ( a1*(s22-s33)**2 + a2*(s33-s11)**2 + a3*(s11-s22)**2 )/6.0 + \
          a4* s23**2      + a5* s31**2      + a6* s12**2

  J30 = ( (b1    +b2    )*s11**3  + (b3    +b4    )*s22**3  + ( b1+b4-b2      +  b1+b4-b3     )*s33**3)/27.0- \
        ( (b1*s22+b2*s33)*s11**2  + (b3*s33+b4*s11)*s22**2  + ((b1+b4-b2)*s11 + (b1+b4-b3)*s22)*s33**2)/9.0 + \
        ( (b1+b4)*s11*s22*s33/9.0 + b11*s12*s23*s31 )*2.0   - \
        ( ( 2.0*b9 *s22 - b8*s33  - (2*b9 -b8)*s11 )*s31**2 + 
          ( 2.0*b10*s33 - b5*s22  - (2*b10-b5)*s11 )*s12**2 +
          ( (b6+b7)*s11 - b6*s22  - b7*s33         )*s23**2
        )/3.0
  f0 = (J20**3 - c*J30**2)**(1.0/6.0)
  k2 = (sigma0/3.0) *18.0 **(1.0/6.0)
  r  = f0/k2 - 1.0
  return r.ravel()

def Cazacu_Barlat2D(sigma0,a1,a2,a3,a6, b1,b2,b3,b4,b5,b10, c,
                    ydata, sigmas):
  '''
    residuum of the Cazacu–Barlat (CZ) yield criterion for plain stress
  '''
  s11 = sigmas[0]; s22 = sigmas[1]; s12 = sigmas[3]

  J20 = ( (a2+a3)*s11**2 + (a1+a3)*s22**2 - 2.0*a3*s11*s22 )/6.0 + a6*s12**2

  J30 = ( (b1     + b2    )*s11**3  + (b3    +b4    )*s22**3 )/27.0- \
        ( (b1*s11 + b4*s22)*s11*s22  )/9.0 + \
        (  b5*s22 + (2*b10-b5)*s11 )*s12**2/3.0
  f0 = (J20**3 - c*J30**2)**(1.0/6.0)
  k2 = (sigma0/3.0) *18.0 **(1.0/6.0)
  r  = f0/k2 - 1.0
  return r.ravel()

def BBC2003(sigma0, a,b,c, d,e,f,g, k, ydata, sigmas):
  '''
    residuum of the BBC2003 yield criterion for plain stress
  '''
  s11 = sigmas[0]; s22 = sigmas[1]; s12 = sigmas[3]
  k2  = 2.0*k

  Gamma =     s11*(d+e)     + s22*(e+f)
  Psi   = ( ( s11*(d-e)/2.0 + s22*(e-f)/2.0 )**2 + (g*s12)**2 )**0.5

  sBar  = ( a*(b*Gamma + c*Psi)**k2 + a*(b*Gamma - c*Psi)**k2 + 
            (1-a)*(2.0*c*Psi)**k2 )**(1.0/k2)
  r = sBar/sigma0 - 1.0
  return r.ravel()

def DruckerBasis(sigma0, C_D, p, ydata, sigmas):
  s11 = sigmas[0]; s22 = sigmas[1]; s33 = sigmas[2]
  s12 = sigmas[3]; s23 = sigmas[4]; s31 = sigmas[5]
  I1  = s11 + s22 + s33
  I2  = s11*s22 + s22*s33 + s11*s33 - s12**2 - s23**2 - s31**2
  I3  = s11*s22*s33 + 2.0*s12*s23*s31 - s12**2*s33 - s23**2*s11 - s31**2*s22
  J2  = I1**2/3.0 - I2
  J3  = I1**3/13.5 - I1*I2/3.0 + I3
  r  = (J2**(3.0*p) - C_D*J3**(2.0*p))*27/(sigma0**6.0) - 1.0
  return r.ravel()

def HosfordBasis(sigma0, F,G,H, a, sigmas, Jac=False, nParas=1):
  '''
    residuum of Hershey yield criterion (eq. 2.43, Y = sigma0)
  '''
  lambdas = principalStresses(sigmas)
  diff23  = abs(lambdas[1,:] - lambdas[2,:])
  diff31  = abs(lambdas[2,:] - lambdas[0,:])
  diff12  = abs(lambdas[0,:] - lambdas[1,:])
  base    = F*diff23**a + G*diff31**a + H*diff12**a; expo = 1.0/a
  left    = base**expo;    right = 2.0**expo*sigma0
  
  if not Jac:
    if nParas == 1: return (left - right).ravel()
    else:           return (left/right - 1.0).ravel()
  else:
    if nParas > 1:
      ln   = lambda x : np.log(x + 1.0e-32)
      dbda = F*ln(diff23)*diff23**a + G*ln(diff31)*diff31**a + H*ln(diff12)*diff12**a
      deda = -expo*expo; dldb = expo*left/base; drda = sigma0*(2.0**expo)*ln(2.0)*deda
      ones = np.ones(np.shape(sigmas)[1]);      jac  = []
    
    if   nParas == 1:     # von Mises
      return ones*(-2.0**0.5)
    elif nParas == 2:     # isotropic Hosford
      j1 = ones*(-2.0**expo)                                              # d[]/dsigma0
      j2 = dldb*dbda + left*ln(base)*deda - drda                          # d[]/da
      for a,b in zip(j1, j2): jac.append([a,b])
      return np.array(jac)
    elif nParas == 5:           # anisotropic Hosford
      j1 = -left/right/sigma0   #ones*(-2.0**expo)                                               # d[]/dsigma0
      j2 = dldb*diff23**a/right; j3 = dldb*diff31**a/right; j4 = dldb*diff12**a/right
      j5 =(dldb*dbda + left*ln(base)*deda)/right + left*(-right**(-2))*drda    # d[]/da
      for a,b,c,d,e in zip(j1, j2,j3,j4,j5): jac.append([a,b,c,d,e])
      return np.array(jac)
      
def Barlat1991Basis(sigma0, a, b, c, f, g, h, order, sigmas):
  '''
    residuum of Barlat 1997 yield criterion
  '''
  cos = np.cos; pi = np.pi; abs = np.abs
  A = a*(sigmas[1] - sigmas[2])
  B = b*(sigmas[2] - sigmas[0])
  C = c*(sigmas[0] - sigmas[1])
  F = f* sigmas[4]
  G = g* sigmas[5]
  H = h* sigmas[3]

  I2 =   (F*F + G*G  +  H*H)/3.0  + ((A-C)**2+(C-B)**2+(B-A)**2)/54.0
  I3 =   (C-B)*(A-C) * (B-A)/54.0 + F*G*H - \
       ( (C-B)*F*F   + (A-C)*G*G  + (B-A)*H*H )/6.0
  theta = np.arccos(I3/I2**1.5)
  Phi = np.sqrt(3.0*I2)* (
        (abs(2.0*cos((2.0*theta + pi)/6.0)))**order +
        (abs(2.0*cos((2.0*theta + pi*3.0)/6.0)))**order +
        (abs(2.0*cos(( 2.0*theta + pi*5.0)/6.0)))**order
                         )**(1.0/order)
#  r   = Phi/2.0**(1.0/order) - sigma0
  r   = Phi/2.0**(1.0/order)/sigma0 - 1.0
  # Phi = (3.0*I2)**(order/2.0) * (
        # (abs(2.0*cos((2.0*theta + pi)/6.0)))    **order +
        # (abs(2.0*cos((2.0*theta + pi*3.0)/6.0)))**order +
        # (abs(2.0*cos((2.0*theta + pi*5.0)/6.0)))**order
                                # )
  # r   = (Phi - 2.0*sigma0**order)**(1.0/order)

  return r.ravel()

fittingCriteria = {
  'tresca'         :{'func' : Tresca,
                     'num'  : 1,'err':np.inf,
                     'name' : 'Tresca',
                     'paras': 'Initial yield stress:',
                     'text' : '\nCoefficient of Tresca criterion:\nsigma0: ',
                     'error': 'The standard deviation error is: '
                    },
  'vonmises'       :{'func' : vonMises,
                     'num'  : 1,'err':np.inf,
                     'name' : 'Huber-Mises-Hencky(von Mises)',
                     'paras': 'Initial yield stress:',
                     'text' : '\nCoefficient of Huber-Mises-Hencky criterion:\nsigma0: ',
                     'error': 'The standard deviation error is: '
                    },
  'hosfordiso'     :{'func' : Hosford,
                     'num'  : 2,'err':np.inf,
                     'name' : 'Gerenal isotropic Hosford',
                     'paras': 'Initial yield stress, a:',
                     'text' : '\nCoefficients of Hosford criterion:\nsigma0, a: ',
                     'error': 'The standard deviation errors are: '
                    },
  'hosfordaniso'   :{'func' : generalHosford,
                     'num'  : 5,'err':np.inf,
                     'name' : 'Gerenal isotropic Hosford',
                     'paras': 'Initial yield stress, F, G, H, a:',
                     'text' : '\nCoefficients of Hosford criterion:\nsigma0, F, G, H, a: ',
                     'error': 'The standard deviation errors are: '
                    },
  'hill1948'       :{'func' : Hill1948,
                     'num'  : 6,'err':np.inf,
                     'name' : 'Hill1948',
                     'paras': 'Normalized [F, G, H, L, M, N]:',
                     'text' : '\nCoefficients of Hill1948 criterion:\n[F, G, H, L, M, N]:'+' '*16,
                     'error': 'The standard deviation errors are: '
                    },
  'drucker'        :{'func' : Drucker,
                     'num'  : 2,'err':np.inf,
                     'name' : 'Drucker',
                     'paras': 'Initial yield stress, C_D:',
                     'text' : '\nCoefficients of Drucker criterion:\nsigma0, C_D: ',
                     'error': 'The standard deviation errors are: '
                    },
  'gdrucker'       :{'func' : generalDrucker,
                     'num'  : 3,'err':np.inf,
                     'name' : 'General Drucker',
                     'paras': 'Initial yield stress, C_D, p:',
                     'text' : '\nCoefficients of Drucker criterion:\nsigma0, C_D, p: ',
                     'error': 'The standard deviation errors are: '
                    },
  'barlat1991iso'  :{'func' : Barlat1991iso,
                     'num'  : 2,'err':np.inf,
                     'name' : 'Barlat1991iso',
                     'paras': 'Initial yield stress, m:',
                     'text' : '\nCoefficients of isotropic Barlat 1991 criterion:\nsigma0, m:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'barlat1991aniso':{'func' : Barlat1991aniso,
                     'num'  : 8,'err':np.inf,
                     'name' : 'Barlat1991aniso',
                     'paras': 'Initial yield stress, a, b, c, f, g, h, m:',
                     'text' : '\nCoefficients of anisotropic Barlat 1991 criterion:\nsigma0, a, b, c, f, g, h, m:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'bbc2003'        :{'func' : BBC2003,
                     'num'  : 9,'err':np.inf,
                     'name' : 'Banabic-Balan-Comsa 2003',
                     'paras': 'Initial yield stress, a, b, c, d, e, f, g, k:',
                     'text' : '\nCoefficients of anisotropic Barlat 1991 criterion:\nsigma0, a, b, c, d, e, f, g, k:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'Cazacu_Barlat2D':{'func' : Cazacu_Barlat2D,
                     'num'  : 12,'err':np.inf,
                     'name' : 'Cazacu Barlat for plain stress',
                     'paras': 'Initial yield stress, a1,a2,a3,a6; b1,b2,b3,b4,b5,b10; c:',
                     'text' : '\nCoefficients of Cazacu Barlat yield criterion for plane stress: \
                               \n Y, a1,a2,a3,a6; b1,b2,b3,b4,b5,b10; c:\n',
                     'error': 'The standard deviation errors are: '
                    },
  'Cazacu_Barlat3D':{'func' : Cazacu_Barlat3D,
                     'num'  : 19,'err':np.inf,
                     'name' : 'Cazacu Barlat',
                     'paras': 'Initial yield stress, a1,a2,a3,a4,a5,a6; b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11; c:',
                     'text' : '\nCoefficients of Cazacu Barlat yield criterion for plane stress: \
                               \n Y, a1,a2,a3,a4,a5,a6; b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11; c\n',
                     'error': 'The standard deviation errors are: '
                    },
  'worst'          :{'err':np.inf},
  'best'           :{'err':np.inf}
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
  def __init__(self,finalStrain,incs,time):
    print('using the random load case generator')
    self.finalStrain = finalStrain
    self.incs = incs
    self.time = time

  def getLoadcase(self,N=0):
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
    criteriaClass = fittingCriteria[nameCriterion]['func'];  criteria = criteriaClass()
    numParas      = fittingCriteria[nameCriterion]['num']
    textParas     = fittingCriteria[nameCriterion]['text'] + formatOutput(numParas)
    textError     = fittingCriteria[nameCriterion]['error']+ formatOutput(numParas,'%-14.8f')+'\n'
    bounds        = fittingCriteria[nameCriterion]['bound']  # Default bounds, no bound
    guess0        = fittingCriteria[nameCriterion]['guess']  # Default initial guess, depends on bounds

    if fitResults == [] : initialguess = guess0
    else                : initialguess = np.array(fitResults[-1])
    weight = get_weight(np.shape(stress)[1])
    ydata  = np.zeros(np.shape(stress)[1])
    try:
      popt, pcov, infodict, errmsg, ierr = \
         leastsqBound (criteria.fun,  initialguess,        args=(ydata,stress),
                       bounds=bounds, full_output=True)
      if ierr not in [1, 2, 3, 4]: raise RuntimeError("Optimal parameters not found: " + errmsg)
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
  me=getLoadcase()
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
  table.data_readArray(['%i_Cauchy'%(i+1) for i in xrange(9)]+[thresholdKey])

  line = 0
  lines = np.shape(table.data)[0]
  yieldStress = np.empty((int(options.yieldValue[2]),6),'d')
  for i,threshold in enumerate(np.linspace(options.yieldValue[0],options.yieldValue[1],options.yieldValue[2])):
    while line < lines:
      if table.data[line,9]>= threshold:
        upper,lower = table.data[line,9],table.data[line-1,9]                                       # values for linear interpolation
        stress = np.array(table.data[line-1,0:9] * (upper-threshold)/(upper-lower) + \
                          table.data[line  ,0:9] * (threshold-lower)/(upper-lower)).reshape(3,3)    # linear interpolation of stress values
        yieldStress[i,0]= stress[0,0]; yieldStress[i,1]=stress[1,1]; yieldStress[i,2]=stress[2,2]
        yieldStress[i,3]=(stress[0,1] + stress[1,0])/2.0     #   0  3  5
        yieldStress[i,4]=(stress[1,2] + stress[2,1])/2.0     #   *  1  4  yieldStress
        yieldStress[i,5]=(stress[2,0] + stress[0,2])/2.0     #   *  *  2
        break
      else:
        line+=1
  
  s.acquire()
  global stressAll
  print('number of yield points of sim %i: %i'%(me,len(yieldStress)))
  print('starting fitting for sim %i from %s'%(me,thread))
  try:
    for i in xrange(int(options.yieldValue[2])):
      stressAll[i]=np.append(yieldStress[i]/unitGPa,stressAll[i])
      myFit.fit(stressAll[i].reshape(len(stressAll[i])//6,6).transpose())
  except Exception as detail:
    print('could not fit for sim %i from %s'%(me,thread))
    print detail
    s.release()
    return
  s.release()

def getLoadcase():
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
parser.set_defaults(min        = 12)
parser.set_defaults(max        = 30)
parser.set_defaults(threads    = 4)
parser.set_defaults(yieldValue = (0.002,0.004,2))
parser.set_defaults(load       = (0.010,100,100.0))
parser.set_defaults(criterion  = 'worst')
parser.set_defaults(fitting    = 'totalshear')
parser.set_defaults(geometry   = '20grains16x16x16')

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

unitGPa = 10.e8
N_simulations=0
fitResults = []
s=threading.Semaphore(1)

stressAll=[np.zeros(0,'d').reshape(0,0) for i in xrange(int(options.yieldValue[2]))]
myLoad = Loadcase(options.load[0],options.load[1],options.load[2])
myFit = Criterion(options.criterion)

threads=[]

for i in range(options.threads):
  threads.append(myThread(i))
  threads[i].start()

for i in range(options.threads):
  threads[i].join()

print 'finished fitting to yield criteria'
