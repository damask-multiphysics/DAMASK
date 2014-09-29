#!/usr/bin/python
# -*- coding: UTF-8 no BOM -*-

import threading,time,os,subprocess,shlex,string,sys,random
import numpy as np
from optparse import OptionParser
from operator import mul
from cStringIO import StringIO
import damask

scriptID = '$Id: yieldSurface.py 3390 2014-08-18 20:09:09Z MPIE\m.diehl $'
scriptName = scriptID.split()[1]
mismatch = None
currentSeedsName = None

def execute(cmd,streamIn=None,dir='./'):
      
  initialPath=os.getcwd()
  os.chdir(dir)
  process = subprocess.Popen(shlex.split(cmd),stdout=subprocess.PIPE,stderr = subprocess.PIPE,stdin=subprocess.PIPE)
  if streamIn != None:
    out,error = process.communicate(streamIn.read())
  else:
    out,error = process.communicate()
  os.chdir(initialPath)

  return out,error

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
    global bestSeedsUpdate
    global bestSeedsVFile
    global nMicrostructures
    global delta
    global points
    global target
    global match
    global baseFile

    s.acquire()
    bestMatch = match
    s.release()
    
    random.seed(options.randomSeed+self.threadID)
    knownSeedsUpdate = bestSeedsUpdate -1.0
    randReset = True
    
    myBestSeedsVFile    = StringIO()
    perturbedSeedsVFile = StringIO()
    perturbedGeomVFile  = StringIO()

    while bestMatch < options.threshold:
# if a better seed file exist, read it into a virtual file
      s.acquire()
      if bestSeedsUpdate > knownSeedsUpdate:
        #print 'Thread %i: '%self.threadID+'Reading new best seeds from', bestSeedsUpdate
        knownSeedsUpdate = bestSeedsUpdate
        bestSeedsVFile.reset()
        myBestSeedsVFile.close()
        myBestSeedsVFile = StringIO()
        i=0
        for line in bestSeedsVFile:
          myBestSeedsVFile.write(line)
      s.release()
# new direction if current one is not good
      if randReset:
        selectedMs = random.randrange(1,nMicrostructures)
        direction = np.array(((random.random()-0.5)*delta[0],
                              (random.random()-0.5)*delta[1],
                              (random.random()-0.5)*delta[2]))
        randReset = False
# set virtual file pointer to the beginning
      perturbedSeedsVFile.close()
      perturbedSeedsVFile = StringIO()
      myBestSeedsVFile.reset()
# perturbe current best fitting seed file  (virtual file lastSeedsVFile)
      perturbedSeedsTable = damask.ASCIItable(myBestSeedsVFile,
                                              perturbedSeedsVFile,
                                              labels=True)
      perturbedSeedsTable.head_read()                                                                                 # read ASCII header info
      perturbedSeedsTable.head_write()
      outputAlive=True
      ms = 1
      while outputAlive and perturbedSeedsTable.data_read():   
        if ms == selectedMs:
          direction+=direction
          newCoords=np.array(tuple(map(float,perturbedSeedsTable.data[0:3]))+direction)
          newCoords=np.where(newCoords>=1.0,newCoords-1.0,newCoords)
          newCoords=np.where(newCoords <1.0,newCoords+1.0,newCoords)
          #print ms,perturbedSeedsTable.data[0:3],'-->',[format(f, '8.6f') for f in newCoords]
          perturbedSeedsTable.data[0:3]=[format(f, '8.6f') for f in newCoords]
        ms+=1
        perturbedSeedsTable.data_write()

# do tesselation with perturbed seed file and evaluate geometry
      perturbedGeomVFile.close()
      perturbedGeomVFile = StringIO()
      perturbedSeedsVFile.reset()
      perturbedGeomVFile.write(execute('geom_fromVoronoiTessellation '+
                     ' -g '+' '.join(map(str, options.grid)),streamIn=perturbedSeedsVFile)[0])
      perturbedGeomVFile.reset()
      perturbedGeomTable = damask.ASCIItable(perturbedGeomVFile,labels=False)
      perturbedGeomTable.head_read()
      for i in perturbedGeomTable.info:
        if i.startswith('microstructures'): myNmicrostructures = int(i.split('\t')[1])
      perturbedGeomTable.data_readArray()
      perturbedGeomTable.output_flush()
      currentData=np.bincount(perturbedGeomTable.data.astype(int).ravel())[1:]/points
      currentError=[]
      currentHist=[]
      for i in xrange(nMicrostructures):
        currentHist.append(np.histogram(currentData,bins=target[i]['bins'])[0])
        currentError.append(np.sqrt(np.square(np.array(target[i]['histogram']-currentHist[i])).sum()))

      s.acquire()
      bestMatch = match
      #print 'Thread %i: finished tessellation (%i microstructures)'%(self.threadID, myNmicrostructures)
      myMatch=0
      for i in xrange(nMicrostructures):
        if currentError[i] > 0.0: break
        myMatch = i+1

      if myNmicrostructures != nMicrostructures:
        print 'Thread %i: Microstructure mismatch (%i microstructures mapped)'%(self.threadID,myNmicrostructures)
        randReset = True
      else:
        for i in xrange(nMicrostructures):

          if currentError[i] > target[i]['error']:
            #print '%i decreased quality to %8.6f'%(i+1,ratio)
            #print 'target',target[i]['histogram']
            #print 'current',currentHist[i]
            randReset = True
            break
          elif currentError[i] < target[i]['error']:
            bestSeedsUpdate = time.time()
            print 'Thread %i: Better match (%i bins, %6.4f --> %6.4f)'%(self.threadID,i+1,target[i]['error'],currentError[i])
            print '          target: ',target[i]['histogram']
            print '          best:   ',currentHist[i]
            currentSeedsName = baseFile+'_'+str(bestSeedsUpdate).replace('.','-')
            perturbedSeedsVFile.reset()
            bestSeedsVFile.close()
            bestSeedsVFile = StringIO()
            sys.stdout.flush()
            with open(currentSeedsName+'.seeds','w') as currentSeedsFile:
              for line in perturbedSeedsVFile:
                currentSeedsFile.write(line)
                bestSeedsVFile.write(line)
            for j in xrange(nMicrostructures):
              target[j]['error'] = currentError[j]
            if myMatch > match:
              print 'Stage %i cleared'%(myMatch)
              match=myMatch
              sys.stdout.flush()
            break
          if i == nMicrostructures-1:
            print 'Thread %i: Continue along trajectory'%(self.threadID)
      
      s.release()


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Monte Carlo simulation to produce seed file that gives same size distribution like given geometry file.

""", version=string.replace(scriptID,'\n','\\n')
)

parser.add_option('-s','--seeds',    dest='seedFile', metavar='string',
                                     help='name of the intial seed file. If not found, a new one is created [%default]')
parser.add_option('-g','--grid',     dest='grid', type='int', nargs=3, metavar='int int int',
                                     help='a,b,c grid of hexahedral box [%default]')
parser.add_option('-t','--threads',  dest='threads', type='int', metavar='int',
                                     help='number of parallel executions [%default]')
parser.add_option('-r', '--rnd',     dest='randomSeed', type='int', metavar='int',
                                     help='seed of random number generator [%default]')
parser.add_option('--target',        dest='target', metavar='string',
                                     help='name of the geom file with target distribution [%default]')
parser.add_option('--tolerance',     dest='threshold', type='int', metavar='int',
                                     help='stopping criterion (bin number) [%default]')
parser.add_option('--scale',         dest='scale',type='float', metavar='float',                                     
                                     help='maximum moving distance of perturbed seed in pixel [%default]')


parser.set_defaults(seedFile    = 'seeds')
parser.set_defaults(grid        = (64,64,64))
parser.set_defaults(threads     = 2)
parser.set_defaults(randomSeed  = None)
parser.set_defaults(target      = 'geom')
parser.set_defaults(threshold   = 20)
parser.set_defaults(scale       = 1.0)


options = parser.parse_args()[0]

if options.randomSeed == None:
  options.randomSeed = int(os.urandom(4).encode('hex'), 16)
print 'random seed', options.randomSeed
delta = (options.scale/options.grid[0],options.scale/options.grid[1],options.scale/options.grid[2])
baseFile=os.path.splitext(os.path.basename(options.seedFile))[0]
points = float(reduce(mul,options.grid))


# ----------- calculate target distribution and bin edges
with open(os.path.splitext(os.path.basename(options.target))[0]+'.geom') as targetGeomFile:
  targetGeomTable = damask.ASCIItable(targetGeomFile,labels=False)
  targetGeomTable.head_read()
  for i in targetGeomTable.info:
    if i.startswith('microstructures'): nMicrostructures = int(i.split()[1])
    if i.startswith('grid'):            targetPoints     = np.array(map(float,i.split()[2:7:2])).prod()

  targetGeomTable.data_readArray()
  targetVolFrac = np.bincount(targetGeomTable.data.astype(int).ravel())[1:nMicrostructures+1]/targetPoints
  target=[]
  for i in xrange(1,nMicrostructures+1):
    targetHist,targetBins = np.histogram(targetVolFrac,bins=i) #bin boundaries
    target.append({'histogram':targetHist,'bins':targetBins})

# ----------- create initial seed file or open existing one
bestSeedsVFile = StringIO()
if os.path.isfile(os.path.splitext(options.seedFile)[0]+'.seeds'):
  with open(os.path.splitext(options.seedFile)[0]+'.seeds') as initialSeedFile:
    for line in initialSeedFile: bestSeedsVFile.write(line)
else:
  bestSeedsVFile.write(execute('seeds_fromRandom'+\
                                ' -g '+' '.join(map(str, options.grid))+\
                                ' -r %i'%options.randomSeed+\
                                ' -N '+str(nMicrostructures))[0])
bestSeedsUpdate = time.time()

# ----------- tessellate initial seed file to get and evaluate geom file
bestSeedsVFile.reset()
initialGeomVFile = StringIO()
initialGeomVFile.write(execute('geom_fromVoronoiTessellation '+
                               ' -g '+' '.join(map(str, options.grid)),bestSeedsVFile)[0])
initialGeomVFile.reset()
initialGeomTable = damask.ASCIItable(initialGeomVFile,labels=False)
initialGeomTable.head_read()
for i in initialGeomTable.info:
  if i.startswith('microstructures'): initialMicrostructures = int(i.split('\t')[1])
if initialMicrostructures != nMicrostructures: print 'error. Microstructure count mismatch'
initialGeomTable.data_readArray()
initialData = np.bincount(initialGeomTable.data.astype(int).ravel())[1:]/points
for i in xrange(nMicrostructures):
  initialHist = np.histogram(initialData,bins=target[i]['bins'])[0]
  target[i]['error']=np.sqrt(np.square(np.array(target[i]['histogram']-initialHist)).sum())
  #print target[i]['histogram']
  #print initialHist
  #print target[i]['error'],'\n-----------------------------------'

match=0
for i in xrange(nMicrostructures):
  if target[i]['error'] > 0.0: break
  match = i+1

print 'Stage %i cleared'%match
initialGeomVFile.close()



# strart mulithreaded monte carlo simulation
threads=[]
s=threading.Semaphore(1)

for i in range(options.threads):
  threads.append(myThread(i))
  threads[i].start()
for i in range(options.threads):
  threads[i].join()
