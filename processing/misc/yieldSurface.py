#!/usr/bin/python
# -*- coding: UTF-8 no BOM -*-

import threading,time,os,subprocess,shlex,string
import numpy as np
from scipy.optimize import curve_fit
from scipy.linalg import svd
from optparse import OptionParser
import damask

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
    eigenvalues = np.linalg.eigvalsh(np.array(sigmas[:,i]).reshape(3,3))
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

def Tresca(sigmas, sigma0):
  '''
    residuum of Tresca yield criterion (eq. 2.26)
  '''
  lambdas = principalStresses(sigmas)
  r = np.amax(np.array([abs(lambdas[2,:]-lambdas[1,:]),\
                        abs(lambdas[1,:]-lambdas[0,:]),\
                        abs(lambdas[0,:]-lambdas[2,:])]),0) - sigma0
  return r.ravel()


def vonMises(sigmas, sigma0):
  '''
    residuum of Huber-Mises-Hencky yield criterion (eq. 2.37)
  '''

  return Hosford(sigmas, sigma0, 2.0)


def Drucker(sigmas, sigma0, C_D):
  '''
    residuum of Drucker yield criterion (eq. 2.41, F = sigma0)
  '''

  return generalDrucker(sigmas, sigma0, C_D, 1.0)


def generalDrucker(sigmas, sigma0, C_D, p):
  '''
    residuum of general Drucker yield criterion (eq. 2.42, F = sigma0)
  '''
  Is = stressInvariants(principalStresses(sigmas))
  r  = (Is[1,:]**(3.0*p)-C_D*Is[2,:]**(2.0*p)) - sigma0
  return r.ravel()


def Hosford(sigmas, sigma0, a):
  '''
    residuum of Hershey yield criterion (eq. 2.43, Y = sigma0)
  '''
  lambdas = principalStresses(sigmas)
  r   = (abs(lambdas[2,:]-lambdas[1,:]))**a\
      + (abs(lambdas[1,:]-lambdas[0,:]))**a\
      + (abs(lambdas[0,:]-lambdas[2,:]))**a\
      -2.0*(abs(sigma0))**a
  return r.ravel()

#more to do
# KarafillisAndBoyce

# ---------------------------------------------------------------------------------------------
# isotropic yield surfaces
# ---------------------------------------------------------------------------------------------

def Hill1948(sigmas, F,G,H,L,M,N):
  '''
    residuum of Hill 1948 quadratic yield criterion (eq. 2.48)
  '''
  r =     F*(sigmas[4]-sigmas[8])**2.0\
    +     G*(sigmas[8]-sigmas[0])**2.0\
    +     H*(sigmas[0]-sigmas[4])**2.0\
    + 2.0*L* sigmas[5]**2.0\
    + 2.0*M* sigmas[2]**2.0\
    + 2.0*N* sigmas[1]**2.0\
    - 1.0
  return r.ravel()/2.0

#more to do
# Hill 1979

# Hill 1990,1993 need special stresses to fit

def generalHosford(sigmas, sigma0, a):
  '''
    residuum of Hershey yield criterion (eq. 2.104, sigma = sigma0)
  '''
  lambdas = principalStresses(sigmas)
  r = np.amax(np.array([F*(abs(lambdas[:,1]-lambdas[:,2]))**a,\
                        G*(abs(lambdas[:,2]-lambdas[:,0]))**a,\
                        H*(abs(lambdas[:,0]-lambdas[:,1]))**a]),1) - sigma0**a
  return r.ravel()


def Barlat1991(sigmas, sigma0, order, \
              a=1.0, b=1.0, c=1.0, f=1.0, g=1.0, h=1.0):
  '''
    residuum of Barlat 1997 yield criterion
  '''
  cos = np.cos; pi = np.pi; abs = np.abs
  A = a*(sigmas[4] - sigmas[8])
  B = b*(sigmas[8] - sigmas[0])
  C = c*(sigmas[0] - sigmas[4])
  F = f*sigmas[5]
  G = g*sigmas[2]
  H = h*sigmas[1]
  I2 = (F*F + G*G + H*H)/3.0 + ((A-C)*(A-C)+(C-B)*(C-B)+(B-A)*(B-A))/54.0
  I3 = (C-B)*(A-C)*(B-A)/54.0 + F*G*H - \
       ((C-B)*F*F+(A-C)*G*G+(B-A)*H*H)/6
  theta = np.arccos(I3/pow(I2,1.5))
  Phi = pow(3.0*I2, order/2.0)* (
    pow(abs(2.0*cos((2.0*theta + pi)/6.0)), order) +
    pow(abs(2.0*cos((2.0*theta + pi*3.0)/6.0)), order) +
    pow(abs(2.0*cos((2.0*theta + pi*5.0)/6.0)), order)
                                 )
  r   = Phi - 2.0*pow(sigma0, order)

  return r.ravel()

def Barlat1991iso(sigmas, sigma0, m):
  '''
    residuum of isotropic Barlat 1991 yield criterion (eq. 2.37)
  '''
  return Barlat1991(sigmas, sigma0, m)

def Barlat1991aniso(sigmas, sigma0, m, a,b,c,f,g,h):
  '''
    residuum of anisotropic Barlat 1991 yield criterion (eq. 2.37)
  '''
  return Barlat1991(sigmas, sigma0, m, a,b,c,f,g,h)


def Barlat1994(sigmas, sigma0, a):
  '''
    residuum of Hershey yield criterion (eq. 2.104, sigma_e = sigma0)
  '''

  return None



fittingCriteria = {
   'Tresca'         :{'fit'  :np.ones(1,'d'),'err':np.inf,
                      'name' :'Tresca',
                               'paras':'Initial yield stress:'},
   'vonMises'       :{'fit'  :np.ones(1,'d'),'err':np.inf,
                      'name' :'Huber-Mises-Hencky(von Mises)',
                      'paras':'Initial yield stress:'},
   'Hosford'        :{'fit'  :np.ones(2,'d'),'err':np.inf,
                      'name' :'Gerenal Hosford',
                      'paras':'Initial yield stress:'},
   'Hill1948'       :{'fit'  :np.ones(6,'d'),'err':np.inf,
                      'name' :'Hill1948',
                      'paras':'Normalized [F, G, H, L, M, N]'},
   'Drucker'        :{'fit'  :np.ones(2,'d'),'err':np.inf,
                      'name' :'Drucker',
                      'paras':'Initial yield stress, C_D:'},
   'Barlat1991iso'  :{'fit'  :np.ones(1,'d'),'err':np.inf,
                      'name' :'Barlat1991iso',
                      'paras':'Initial yield stress, m:'},
   'Barlat1991aniso':{'fit'  :np.ones(7,'d'),'err':np.inf,
                      'name' :'Barlat1991aniso',
                      'paras':'Initial yield stress, m, a, b, c, f, g, h:'},
   'worst'          :{'err':np.inf},
   'best'           :{'err':np.inf}
                   }


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
    if self.name.lower() == 'tresca':
      funResidum = Tresca
      text = '\nCoefficient of Tresca criterion:\nsigma0: '+formatOutput(1)
      error='The standard deviation error is: '+formatOutput(1,'%-14.8f')+'\n'
    elif self.name.lower() == 'vonmises':
      funResidum = vonMises
      text = '\nCoefficient of Huber-Mises-Hencky criterion:\nsigma0: '+formatOutput(1)
      error='The standard deviation error is: '+formatOutput(1,'%-14.8f')+'\n'
    elif self.name.lower() == 'hosford':
      funResidum = Hosford
      text = '\nCoefficient of general Hosford criterion:\nsigma0, a: '+formatOutput(2)
      error='The standard deviation error is: '+formatOutput(2,'%-14.8f')+'\n'
    elif self.name.lower() == 'drucker':
      funResidum = Drucker
      text = '\nCoefficient of Drucker criterion:\nsigma0, C_D: '+formatOutput(2)
      error='The standard deviation errors are: '+formatOutput(2,'%-14.8f')+'\n'
    elif self.name.lower() == 'hill1948':
      funResidum = Hill1948
      text = '\nCoefficient of Hill1948 criterion:\n[F, G, H, L, M, N]:'+' '*16+formatOutput(6)
      error='The standard deviation errors are: '+formatOutput(6,'%-14.8f')+'\n'
    elif self.name.lower() == 'barlat1991iso':
      funResidum = Barlat1991iso
      text = '\nCoefficient of isotropic Barlat 1991 criterion:\nsigma0, m:\n'+formatOutput(2)
      error='The standard deviation errors are: '+formatOutput(1,'%-14.8f')+'\n'
    elif self.name.lower() == 'barlat1991aniso':
      funResidum = Barlat1991aniso
      text = '\nCoefficient of anisotropic Barlat 1991 criterion:\nsigma0, \m, a, b, c, f, g, h:\n' \
             +formatOutput(8)
      error='The standard deviation errors are: '+formatOutput(7,'%-14.8f')
    if fitResults == []:
      initialguess = fittingCriteria[funResidum.__name__]['fit']
    else:
      initialguess = np.array(fitResults[-1])
    weight = get_weight(np.shape(stress)[1])
    try:
      popt, pcov = \
        curve_fit(funResidum, stress, np.zeros(np.shape(stress)[1]),
                  initialguess, weight)
      perr = np.sqrt(np.diag(pcov))
      fitResults.append(popt.tolist())
      print (text%array2tuple(popt))
      print (error%array2tuple(perr))
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
  yieldStress=[None for i in xrange(int(options.yieldValue[2]))]
  for i,threshold in enumerate(np.linspace(options.yieldValue[0],options.yieldValue[1],options.yieldValue[2])):
    while line < lines:
      if table.data[line,9]>= threshold:
        upper,lower = table.data[line,9],table.data[line-1,9]                                       # values for linear interpolation
        yieldStress[i] = table.data[line-1,0:9] * (upper-threshold)/(upper-lower) \
                       + table.data[line  ,0:9] * (threshold-lower)/(upper-lower)                   # linear interpolation of stress values
        yieldStress[i][1] = (yieldStress[i][1] + yieldStress[i][3])/2.0
        yieldStress[i][2] = (yieldStress[i][2] + yieldStress[i][6])/2.0
        yieldStress[i][5] = (yieldStress[i][5] + yieldStress[i][7])/2.0
        yieldStress[i][3] = yieldStress[i][1]
        yieldStress[i][6] = yieldStress[i][2]
        yieldStress[i][7] = yieldStress[i][5]
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
      myFit.fit(stressAll[i].reshape(len(stressAll[i])//9,9).transpose())
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
