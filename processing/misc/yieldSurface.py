#!/usr/bin/python
# -*- coding: UTF-8 no BOM -*-

import numpy as np
from scipy.optimize import curve_fit
from scipy.linalg import svd
import threading,time,os,subprocess,shlex
import damask
geomName='20grains16x16x16'
popt1=np.ones(6,'d')

def execute(cmd,dir='./'):
      
  initialPath=os.getcwd()
  out = ''
  line = True
  process = subprocess.Popen(shlex.split(cmd),stdout=subprocess.PIPE,stderr = subprocess.STDOUT)
  while line:
    line = process.stdout.readline()
    out += line
  os.chdir(initialPath)

def asFullTensor(voigt):
  return np.array([[voigt[0],voigt[3],voigt[5]],\
          [voigt[3],voigt[1],voigt[4]],\
          [voigt[5],voigt[4],voigt[2]]])

def Hill48(x, F,G,H,L,M,N):
  a= F*(x[1]-x[2])**2 + G*(x[2]-x[0])**2 + H*(x[0]-x[1])** + \
         2*L*x[4]**2 + 2*M*x[5]**2 + 2*N*x[3]**2 -1.
  return a.ravel()

def vonMises(x, S_y):
  sv=np.zeros(0,'d')
  for i in xrange(np.shape(x)[0]):
    U, l, Vh = svd(np.array(x[i,:]).reshape(3,3))
    sv = np.append(sv,l)
  sv = sv.reshape(np.shape(x)[0],3)
  ooo= (sv[:,2]-sv[:,1])**2+(sv[:,1]-sv[:,0])**2+(sv[:,0]-sv[:,2])**2-2*S_y**2
  return ooo

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
    for i in main[:2]: # fill 2 out of 3 main entries
      defgrad[i]=1.+values[i]
      stress[i]='*'
    for off in [[1,3,0],[2,6,0],[5,7,0]]: # fill 3 off-diagonal pairs of defgrad (1 or 2 entries)
      off=np.array(off)
      np.random.shuffle(off)
      for i in off[0:2]:
        if i != 0: 
          defgrad[i]=values[i]
          stress[i]='*'

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
    try:
      spopt1, pcov = curve_fit(Hill48, stress, np.zeros(np.shape(stress)[1]),p0=popt1)
      popt = popt1
      print 'Hill 48', popt
      print Hill48(stress,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])
    except Exception as detail:
      print detail
      pass
    #popt, pcov = curve_fit(vonMises, stress, np.zeros(np.shape(stress)[0]))
    #print 'von Mises', popt


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
  if not os.path.isfile('%s.load'%me):
    print('generating loadcase for sim %s from %s'%(me,thread))
    f=open('%s.load'%me,'w')
    f.write(myLoad.getNext(me))
    f.close()
    s.release()
  else: s.release()
  
  s.acquire()
  if not os.path.isfile('%s_%s.spectralOut'%(geomName,me)):
    print('starting simulation %i from %s'%(me,thread))
    s.release()
    execute('DAMASK_spectral -l %i -g %i'%(me,geomName))
  else: s.release()
   

  s.acquire()
  if not os.path.isfile('./postProc/%s_%s.txt'%(geomName,me)):
    print('starting post processing for sim %i from %s'%(me,thread))
    s.release()
    execute('postResults --cr f,p %i_%i.spectralOut'%(geomName,me))
    execute('addCauchy ./postProc/%i_%i.txt'%(geomName,me))
    execute('addStrainTensors -l -v ./postProc/%i_%i.txt'%(geomName,me))
    execute('addMises -s Cauchy -e ln(V) ./postProc/%i_%i.txt'%(geomName,me))
  else: s.release()

  s.acquire()
  print('reading values for sim %i from %s'%(me,thread))
  s.release()

  refFile = open('./postProc/20grains16x16x16_%i.txt'%me)
  refFile = open('./postProc/20grains16x16x16_1.txt')
  table = damask.ASCIItable(refFile)
  table.head_read()
  for l in ['Mises(ln(V))','1_Cauchy']:
    if l not in table.labels: print '%s not found'%l
  while table.data_read():
    if float(table.data[table.labels.index('Mises(ln(V))')]) > 0.002:
      yieldStress = np.array(table.data[table.labels.index('1_Cauchy'):table.labels.index('9_Cauchy')+1],'d')
  
  s.acquire()
  print('starting fitting for sim %i from %s'%(me,thread))
  global stressAll
  stressAll=np.append(yieldStress,stressAll)
  for i in range(np.shape(stressAll.reshape(len(stressAll)//9,9).transpose())[1]):
    print i, stressAll.reshape(len(stressAll)//9,9).transpose()[:,i]
  myFit.fit(stressAll.reshape(len(stressAll)//9,9).transpose())
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
maxN_simulations=14
N_simulations=0
s=threading.Semaphore(1)
scale = 0.02
incs = 10
duration = 10
stressAll=np.zeros(0,'d').reshape(0,0)

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
