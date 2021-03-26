#!/usr/bin/env python3

import threading
import time
import os
import sys
import random
from optparse import OptionParser
from io import StringIO

import numpy as np

import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

mismatch = None
currentSeedsName = None

#---------------------------------------------------------------------------------------------------
class myThread (threading.Thread):
  """Perturb seed in seed file, performes Voronoi tessellation, evaluates, and updates best match."""

  def __init__(self, threadID):
    """Threading class with thread ID."""
    threading.Thread.__init__(self)
    self.threadID = threadID

  def run(self):
    global bestSeedsUpdate
    global bestSeedsVFile
    global nMaterials
    global delta
    global points
    global target
    global match
    global baseFile
    global maxSeeds

    s.acquire()
    bestMatch = match
    s.release()

    random.seed(options.randomSeed+self.threadID)                                                   # initializes to given seeds
    knownSeedsUpdate = bestSeedsUpdate -1.0                                                         # trigger update of local best seeds
    randReset = True                                                                                # aquire new direction

    myBestSeedsVFile    = StringIO()                                                                # store local copy of best seeds file
    perturbedSeedsVFile = StringIO()                                                                # perturbed best seeds file

#--- still not matching desired bin class ----------------------------------------------------------
    while bestMatch < options.threshold:
      s.acquire()                                                                                   # ensure only one thread acces global data
      if bestSeedsUpdate > knownSeedsUpdate:                                                        # write best fit to virtual file
        knownSeedsUpdate = bestSeedsUpdate
        bestSeedsVFile.seek(0)
        myBestSeedsVFile.close()
        myBestSeedsVFile = StringIO()
        i=0
        myBestSeedsVFile.writelines(bestSeedsVFile.readlines())
      s.release()

      if randReset:                                                                                 # new direction because current one led to worse fit

        randReset = False

        NmoveGrains = random.randrange(1,maxSeeds)
        selectedMs = []
        direction = []
        for i in range(NmoveGrains):
          selectedMs.append(random.randrange(1,nMaterials))

          direction.append((np.random.random()-0.5)*delta)

      perturbedSeedsVFile.close()                                                                   # reset virtual file
      perturbedSeedsVFile = StringIO()
      myBestSeedsVFile.seek(0)

      perturbedSeedsTable = damask.Table.load(myBestSeedsVFile)
      coords = perturbedSeedsTable.get('pos')
      i = 0
      for ms,coord in enumerate(coords):
        if ms in selectedMs:
          newCoords=coord+direction[i]
          newCoords=np.where(newCoords>=1.0,newCoords-1.0,newCoords)                                # ensure that the seeds remain in the box
          newCoords=np.where(newCoords <0.0,newCoords+1.0,newCoords)
          coords[i]=newCoords
          direction[i]*=2.
          i+= 1
      perturbedSeedsTable.set('pos',coords).save(perturbedSeedsVFile,legacy=True)

#--- do tesselation with perturbed seed file ------------------------------------------------------
      perturbedGeom = damask.Grid.from_Voronoi_tessellation(options.grid,np.ones(3),coords)


#--- evaluate current seeds file ------------------------------------------------------------------
      myNmaterials = len(np.unique(perturbedGeom.material))
      currentData = np.bincount(perturbedGeom.material.ravel())[1:]/points
      currentError=[]
      currentHist=[]
      for i in range(nMaterials):                                                             # calculate the deviation in all bins per histogram
        currentHist.append(np.histogram(currentData,bins=target[i]['bins'])[0])
        currentError.append(np.sqrt(np.square(np.array(target[i]['histogram']-currentHist[i])).sum()))

# as long as not all grains are within the range of the target, use the deviation to left and right as error
      if currentError[0]>0.0:
        currentError[0] *=((target[0]['bins'][0]-np.min(currentData))**2.0+
                           (target[0]['bins'][1]-np.max(currentData))**2.0)**0.5                    # norm of deviations by number of usual bin deviation
      s.acquire()                                                                                   # do the evaluation serially
      bestMatch = match
#--- count bin classes with no mismatch ----------------------------------------------------------------------
      myMatch=0
      for i in range(nMaterials):
        if currentError[i] > 0.0: break
        myMatch = i+1

      if myNmaterials == nMaterials:
        for i in range(min(nMaterials,myMatch+options.bins)):
          if currentError[i] > target[i]['error']:                                                  # worse fitting, next try
            randReset = True
            break
          elif currentError[i] < target[i]['error']:                                                # better fit
            bestSeedsUpdate = time.time()                                                           # save time of better fit
            damask.util.croak('Thread {:d}: Better match ({:d} bins, {:6.4f} --> {:6.4f})'\
                                          .format(self.threadID,i+1,target[i]['error'],currentError[i]))
            damask.util.croak('          target: '+np.array_str(target[i]['histogram']))
            damask.util.croak('          best:   '+np.array_str(currentHist[i]))
            currentSeedsName = baseFile+'_'+str(bestSeedsUpdate).replace('.','-')                   # name of new seed file (use time as unique identifier)
            perturbedSeedsVFile.seek(0)
            bestSeedsVFile.close()
            bestSeedsVFile = StringIO()
            sys.stdout.flush()
            with open(currentSeedsName+'.seeds','w') as currentSeedsFile:                           # write to new file
              for line in perturbedSeedsVFile:
                currentSeedsFile.write(line)
                bestSeedsVFile.write(line)
            for j in range(nMaterials):                                                       # save new errors for all bins
              target[j]['error'] = currentError[j]
            if myMatch > match:                                                                     # one or more new bins have no deviation
              damask.util.croak( 'Stage {:d}  cleared'.format(myMatch))
              match=myMatch
              sys.stdout.flush()
            break
          if i == min(nMaterials,myMatch+options.bins)-1:                                     # same quality as before: take it to keep on moving
            bestSeedsUpdate = time.time()
            perturbedSeedsVFile.seek(0)
            bestSeedsVFile.close()
            bestSeedsVFile = StringIO()
            bestSeedsVFile.writelines(perturbedSeedsVFile.readlines())
            for j in range(nMaterials):
              target[j]['error'] = currentError[j]
            randReset = True
      else:                                                                                         #--- not all grains are tessellated
        damask.util.croak('Thread {:d}: Material mismatch ({:d} material indices mapped)'\
                                                         .format(self.threadID,myNmaterials))
        randReset = True


      s.release()


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(usage='%prog options [file[s]]', description = """
Monte Carlo simulation to produce seed file that gives same size distribution like given geometry file.

""", version = scriptID)

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
parser.add_option('--bins',          dest='bins', type='int', metavar='int',
                                     help='bins to sort beyond current best fit [%default]')
parser.add_option('--maxseeds',      dest='maxseeds', type='int', metavar='int',
                                     help='maximum number of seeds to move simulateneously [number of seeds]')

parser.set_defaults(seedFile    = 'seeds',
                    grid        = (64,64,64),
                    threads     = 2,
                    randomSeed  = None,
                    target      = 'geom',
                    threshold   = 20,
                    bins        = 15,
                    scale       = 1.0,
                    maxseeds    = 0)

options = parser.parse_args()[0]

damask.util.report(scriptName,options.seedFile)

if options.randomSeed is None:
  options.randomSeed = int(os.urandom(4).hex(),16)
damask.util.croak(options.randomSeed)
delta = options.scale/np.array(options.grid)
baseFile = os.path.splitext(os.path.basename(options.seedFile))[0]
points = np.array(options.grid).prod().astype('float')

# ----------- calculate target distribution and bin edges
targetGeom = damask.Grid.load_ASCII(os.path.splitext(os.path.basename(options.target))[0]+'.geom')
nMaterials = len(np.unique(targetGeom.material))
targetVolFrac = np.bincount(targetGeom.material.flatten())/targetGeom.cells.prod().astype(np.float)
target = []
for i in range(1,nMaterials+1):
  targetHist,targetBins = np.histogram(targetVolFrac,bins=i) #bin boundaries
  target.append({'histogram':targetHist,'bins':targetBins})

# ----------- create initial seed file or open existing one
bestSeedsVFile = StringIO()
if os.path.isfile(os.path.splitext(options.seedFile)[0]+'.seeds'):
  initial_seeds = damask.Table.load(os.path.splitext(options.seedFile)[0]+'.seeds').get('pos')
else:
  initial_seeds = damask.seeds.from_random(np.ones(3),nMaterials,options.grid,options.randomSeed)

bestSeedsUpdate = time.time()

# ----------- tessellate initial seed file to get and evaluate geom file
bestSeedsVFile.seek(0)
initialGeom = damask.Grid.from_Voronoi_tessellation(options.grid,np.ones(3),initial_seeds)

if len(np.unique(targetGeom.material)) != nMaterials:
  damask.util.croak('error. Material count mismatch')

initialData = np.bincount(initialGeom.material.flatten())/points
for i in range(nMaterials):
  initialHist = np.histogram(initialData,bins=target[i]['bins'])[0]
  target[i]['error']=np.sqrt(np.square(np.array(target[i]['histogram']-initialHist)).sum())

# as long as not all grain sizes are within the range, the error is the deviation to left and right
if target[0]['error'] > 0.0:
  target[0]['error'] *=((target[0]['bins'][0]-np.min(initialData))**2.0+
                        (target[0]['bins'][1]-np.max(initialData))**2.0)**0.5
match=0
for i in range(nMaterials):
  if target[i]['error'] > 0.0: break
  match = i+1


if options.maxseeds < 1:
  maxSeeds = len(np.unique(initialGeom.material))
else:
  maxSeeds = options.maxseeds

if match >0: damask.util.croak('Stage {:d} cleared'.format(match))
sys.stdout.flush()

# start mulithreaded monte carlo simulation
threads = []
s = threading.Semaphore(1)

for i in range(options.threads):
  threads.append(myThread(i))
  threads[i].start()
for i in range(options.threads):
  threads[i].join()
