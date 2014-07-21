#!/usr/bin/python
# -*- coding: UTF-8 no BOM -*-

import numpy as np
from scipy.optimize import curve_fit
from scipy.linalg import svd
import threading,time,os,subprocess,shlex,string
import damask
from optparse import OptionParser

scriptID='aa'
popt1=[np.ones(1,'d'),np.ones(6,'d')]

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
  a = F*(x[1]-x[2])**2 + G*(x[2]-x[0])**2 + H*(x[0]-x[1])** + \
         2*L*x[4]**2 + 2*M*x[5]**2 + 2*N*x[3]**2 -1.
  return a.ravel()

def vonMises(x, S_y):
  sv=np.zeros(0,'d')
  for i in xrange(np.shape(x)[1]):
    U, l, Vh = svd(np.array(x[:,i]).reshape(3,3))
    sv = np.append(sv,l)
  sv = sv.reshape(np.shape(x)[1],3)
  ooo = (sv[:,2]-sv[:,1])**2+(sv[:,1]-sv[:,0])**2+(sv[:,0]-sv[:,2])**2-2*S_y**2
  return ooo.ravel()

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
          ' incs %s'%self.incs+\
          ' time %s'%self.time

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

  def fit(self,stress):
    global popt1
    try:  
      popt1[0], pcov = curve_fit(vonMises, stress, np.zeros(np.shape(stress)[1]),p0=popt1[0])
      print 'Mises', popt1[0], pcov
    except Exception as detail:
      print detail
      pass
    try:  
      popt1[1], pcov = curve_fit(Hill48, stress, np.zeros(np.shape(stress)[1]),p0=popt1[1])
      print 'Hill48', popt1[1], pcov
    except Exception as detail:
      print detail
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
  
  global geomName
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
  if not os.path.isfile('%s_%i.spectralOut'%(geomName,me)):
    print('starting simulation %s from %s'%(me,thread))
    s.release()
    execute('DAMASK_spectral -g %s -l %i'%(geomName,me))
  else: s.release()
   
  s.acquire()
  if not os.path.isfile('./postProc/%s_%i.txt'%(geomName,me)):
    print('starting post processing for sim %i from %s'%(me,thread))
    s.release()
    execute('postResults --cr f,p %s_%i.spectralOut'%(geomName,me))
    execute('addCauchy ./postProc/%s_%i.txt'%(geomName,me))
    execute('addStrainTensors -l -v ./postProc/%s_%i.txt'%(geomName,me))
    execute('addMises -s Cauchy -e ln(V) ./postProc/%s_%i.txt'%(geomName,me))
  else: s.release()

  s.acquire()
  print('reading values for sim %i from %s'%(me,thread))
  s.release()

  refFile = open('./postProc/%s_%i.txt'%(geomName,me))
  table = damask.ASCIItable(refFile)
  table.head_read()
  for l in ['Mises(ln(V))','1_Cauchy']:
    if l not in table.labels: print '%s not found'%l
  while table.data_read():
    if float(table.data[table.labels.index('Mises(ln(V))')]) > 0.002:
      yieldStress = np.array(table.data[table.labels.index('1_Cauchy'):table.labels.index('9_Cauchy')+1],'d')/10.e8
  
  s.acquire()
  print('starting fitting for sim %i from %s'%(me,thread))
  global stressAll
  stressAll=np.append(yieldStress,stressAll)
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

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Performs calculations with various loads on given geometry file and fits yield surface.

""", version=string.replace(scriptID,'\n','\\n')
)

parser.add_option('-l','--load' ,    dest='load', type='float', nargs=3, \
                                     help='load: final strain; increments; time', metavar='float int float')
parser.add_option('-g','--geometry', dest='geometry', type='string', \
                                     help='name of the geometry file', metavar='string')
#parser.add_option('-c','--criterion',dest='formulas', action='extend', type='string', \
#                                     help='(list of) formulas corresponding to labels', metavar='<LIST>')
parser.set_defaults(load    = [0.008,80,80.0])
parser.set_defaults(geometry ='20grains16x16x16')

options = parser.parse_args()[0]

geomName =options.geometry
minN_simulations=20
maxN_simulations=40
N_simulations=0
s=threading.Semaphore(1)
scale = 0.02
stressAll=np.zeros(0,'d').reshape(0,0)

myLoad = Loadcase(options.load[0],options.load[1],options.load[2])
myFit = Criterion('vonmises')

N_threads=4
t=[]

for i in range(N_threads):
  t.append(myThread(i))
  t[i].start()

for i in range(N_threads):
  t[i].join()
print "Exiting Main Thread"
