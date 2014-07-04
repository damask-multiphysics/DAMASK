#!/usr/bin/python
# -*- coding: UTF-8 no BOM -*-

import numpy as np
from scipy.optimize import curve_fit
from scipy.linalg import svd
import threading
import time

def asFullTensor(voigt):
  return np.array([[voigt[0],voigt[3],voigt[5]],\
                   [voigt[3],voigt[1],voigt[4]],\
                   [voigt[5],voigt[4],voigt[2]]])

def Hill48(x, F,G,H,L,M,N):
  a= F*(x[1]-x[2])**2 + G*(x[2]-x[0])**2 + H*(x[0]-x[1])** + \
         2*L*x[4]**2 + 2*M*x[5]**2 + 2*N*x[3]**2 -1.
  return a.ravel()

def vonMises(x, S_y):
  p = svd(asFullTensor(x)) 
  return (s[2]-s[1])**2+(s[1]-s[0])**2+(s[0]-s[2])**2-2*S_y**2

#---------------------------------------------------------------------------------------------------
class Loadcase():
#---------------------------------------------------------------------------------------------------
  '''
     Class for generating load cases for the spectral solver
  '''

# ------------------------------------------------------------------
  def __init__(self):
    print('using the random load case generator')

  def getNext(self,N=0):
    defgrad=['*']*9
    stress =[0]*9
    values=(np.random.random_sample(9)-.5)*scale*2

    main=np.array([0,4,8])
    np.random.shuffle(main)
    for i in main[:2]:
      defgrad[i]=1.+values[i]
      stress[i]='*'
    for off in [[1,3,0],[2,6,0],[5,7,0]]:
      off=np.array(off)
      np.random.shuffle(off)
      if off[0] != 0: 
        defgrad[off[0]]=values[off[0]]
        stress[off[0]]='*'
    return 'f '+' '.join(str(c) for c in defgrad)+\
          ' p '+' '.join(str(c) for c in stress)+\
          ' incs %s'%incs+\
          ' time %s'%duration

#---------------------------------------------------------------------------------------------------
class Criterion(object):
#---------------------------------------------------------------------------------------------------
  '''
     Fitting to certain criterion
  '''
  def __init__(self,name):
    self.name = name.lower()
    if self.name not in ['hill48','vonmises']: print('Mist')
    print('using the %s criterion'%self.name)
    self.popt = 0.0

  def fit(self,stress):
    if self.name == 'hill48':
      try:
        self.popt, pcov = curve_fit(Hill48, stress, np.zeros(np.shape(stress)[1]))
        print self.popt
      except:
        pass
    elif self.name == 'vonmises':
      try:
        self.popt, pcov = curve_fit(vonMises, stress.transpose(), np.shape(stress)[1])
      except:
        pass 


#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
  '''
     Runner class
  '''
class myThread (threading.Thread):
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
  print('starting sim %i from %s'%(me,thread))
  f=open('%s.load'%me,'w')
  f.write(myLoad.getNext(me))
  f.close()
  #dummy
  print('doing postprocessing sim %i from %s'%(me,thread))
  voigt = np.random.random_sample(6)*90e6
  global stressAll
  stressAll=np.append(stressAll,asFullTensor(voigt).reshape(9))
  stressAll=stressAll.reshape(9,len(stressAll)//9)
  myFit.fit(stressAll)
  s.release()
  time.sleep(delay)
  s.acquire()
  s.release()

def getLoadcase():
  global N_simulations
  N_simulations+=1
  return N_simulations

def converged():
  global N_simulations
  global maxN_simulations
  if N_simulations < maxN_simulations:
    return False
  else:
    return True

# main

minN_simulations=20
maxN_simulations=10
N_simulations=0
s=threading.Semaphore(1)
scale = 0.02
incs = 10
duration = 10
stressAll=np.zeros(0,'d')

myLoad = Loadcase()
myFit = Criterion('Hill48')

N_threads=3
t=[]

for i in range(N_threads):
  t.append(myThread(i))
  t[i].start()

for i in range(N_threads):
  t[i].join()

print "Exiting Main Thread"
