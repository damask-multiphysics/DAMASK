#!/usr/bin/env python

from optparse import OptionParser
import damask
import os,sys,math,re,random,string

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = scriptID.split()[1]

random.seed()

# --- helper functions ---

def binAsBins(bin,intervals):
  """ explode compound bin into 3D bins list """
  bins = [0]*3
  bins[0] = (bin//(intervals[1] * intervals[2])) % intervals[0]
  bins[1] = (bin//intervals[2]) % intervals[1]
  bins[2] = bin % intervals[2]
  return bins
  
def binsAsBin(bins,intervals):
  """ implode 3D bins into compound bin """
  return (bins[0]*intervals[1] + bins[1])*intervals[2] + bins[2]

def EulersAsBins(Eulers,intervals,deltas,center):
  """ return list of Eulers translated into 3D bins list """
  return [\
  int((euler+(0.5-center)*delta)//delta)%interval \
  for euler,delta,interval in zip(Eulers,deltas,intervals) \
  ]

def binAsEulers(bin,intervals,deltas,center):
  """ compound bin number translated into list of Eulers """
  Eulers = [0.0]*3
  Eulers[2] = (bin%intervals[2] + center)*deltas[2]
  Eulers[1] = (bin//intervals[2]%intervals[1] + center)*deltas[1]
  Eulers[0] = (bin//(intervals[2]*intervals[1]) + center)*deltas[0]
  return Eulers

def directInvRepetitions(probability,scale):
  """ calculate number of samples drawn by direct inversion """
  nDirectInv = 0
  for bin in range(len(probability)): # loop over bins
    nDirectInv += int(round(probability[bin]*scale)) # calc repetition
  return nDirectInv


# --- sampling methods ---

# ----- efficient algorithm ---------

def directInversion (ODF,nSamples):
  """ ODF contains 'dV_V' (normalized to 1), 'center', 'intervals', 'limits' (in radians) """
  
  nBins = ODF['intervals'][0]*ODF['intervals'][1]*ODF['intervals'][2]
  deltas = [limit/intervals for limit,intervals in zip(ODF['limits'],ODF['intervals'])]
  
  # calculate repetitions of each orientation
  if re.search(r'hybrid',sys.argv[0],re.IGNORECASE):    # my script's name contains "hybrid"
    nOptSamples = max(ODF['nNonZero'],nSamples)         # random subsampling if too little samples requested
  else:                                                 # blunt integer approximation
    nOptSamples = nSamples

  nInvSamples = 0
  repetition = [None]*nBins
  probabilityScale = nOptSamples # guess

  scaleLower = 0.0
  nInvSamplesLower = 0
  scaleUpper = float(nOptSamples)
  incFactor = 1.0
  nIter = 0
  nInvSamplesUpper = directInvRepetitions(ODF['dV_V'],scaleUpper)
  while (\
      (scaleUpper-scaleLower > scaleUpper*1e-15 or nInvSamplesUpper < nOptSamples) and \
      nInvSamplesUpper != nOptSamples \
      ): # closer match required?
    if nInvSamplesUpper < nOptSamples:
      scaleLower,scaleUpper = scaleUpper,scaleUpper+incFactor*(scaleUpper-scaleLower)/2.0
      incFactor *= 2.0
      nInvSamplesLower,nInvSamplesUpper = nInvSamplesUpper,directInvRepetitions(ODF['dV_V'],scaleUpper)
    else:
      scaleUpper = (scaleLower+scaleUpper)/2.0
      incFactor = 1.0
      nInvSamplesUpper = directInvRepetitions(ODF['dV_V'],scaleUpper)
    nIter += 1
    print '%i:(%12.11f,%12.11f) %i <= %i <= %i'%(nIter,scaleLower,scaleUpper,nInvSamplesLower,nOptSamples,nInvSamplesUpper)
  nInvSamples = nInvSamplesUpper
  scale = scaleUpper
  print 'created set of',nInvSamples,'samples (',float(nInvSamples)/nOptSamples-1.0,') with scaling',scale,'delivering',nSamples
  repetition = [None]*nBins # preallocate and clear
    
  for bin in range(nBins): # loop over bins
    repetition[bin] = int(round(ODF['dV_V'][bin]*scale)) # calc repetition

  # build set
  set = [None]*nInvSamples
  i = 0
  for bin in range(nBins):
    set[i:i+repetition[bin]] = [bin]*repetition[bin] # fill set with bin, i.e. orientation
    i += repetition[bin] # advance set counter
  
  orientations = [None]*nSamples
  reconstructedODF = [0.0]*nBins
  unitInc = 1.0/nSamples
  for j in range(nSamples):
    if (j == nInvSamples-1): ex = j
    else: ex = int(round(random.uniform(j+0.5,nInvSamples-0.5)))
    bin = set[ex]
    bins = binAsBins(bin,ODF['intervals'])
    Eulers = binAsEulers(bin,ODF['intervals'],deltas,ODF['center'])
    orientations[j] = '%g\t%g\t%g' %( math.degrees(Eulers[0]),math.degrees(Eulers[1]),math.degrees(Eulers[2]) )
    reconstructedODF[bin] += unitInc
    set[ex] = set[j] # exchange orientations
  
  return orientations, reconstructedODF


# ----- trial and error algorithms ---------

def MonteCarloEulers (ODF,nSamples):
  """ ODF contains 'dV_V' (normalized to 1), 'center', 'intervals', 'limits' (in radians) """
  
  countMC = 0
  maxdV_V = max(ODF['dV_V'])
  nBins = ODF['intervals'][0]*ODF['intervals'][1]*ODF['intervals'][2]
  deltas = [limit/intervals for limit,intervals in zip(ODF['limits'],ODF['intervals'])]
  orientations = [None]*nSamples
  reconstructedODF = [0.0]*nBins
  unitInc = 1.0/nSamples
  
  for j in range(nSamples):
    MC = maxdV_V*2.0
    bin = 0
    while MC > ODF['dV_V'][bin]:
      countMC += 1
      MC = maxdV_V*random.random()
    Eulers = [limit*random.random() for limit in ODF['limits']]
    bins = EulersAsBins(Eulers,ODF['intervals'],deltas,ODF['center'])
    bin = binsAsBin(bins,ODF['intervals'])
    orientations[j] = '%g\t%g\t%g' %( math.degrees(Eulers[0]),math.degrees(Eulers[1]),math.degrees(Eulers[2]) )
    reconstructedODF[bin] += unitInc

  return orientations, reconstructedODF, countMC


def MonteCarloBins (ODF,nSamples):
  """ ODF contains 'dV_V' (normalized to 1), 'center', 'intervals', 'limits' (in radians) """
  
  countMC = 0
  maxdV_V = max(ODF['dV_V'])
  nBins = ODF['intervals'][0]*ODF['intervals'][1]*ODF['intervals'][2]
  deltas = [limit/intervals for limit,intervals in zip(ODF['limits'],ODF['intervals'])]
  orientations = [None]*nSamples
  reconstructedODF = [0.0]*nBins
  unitInc = 1.0/nSamples
  
  for j in range(nSamples):
    MC = maxdV_V*2.0
    bin = 0
    while MC > ODF['dV_V'][bin]:
      countMC += 1
      MC  = maxdV_V*random.random()
      bin = int(nBins * random.random())
    Eulers = binAsEulers(bin,ODF['intervals'],deltas,ODF['center'])
    orientations[j] = '%g\t%g\t%g' %( math.degrees(Eulers[0]),math.degrees(Eulers[1]),math.degrees(Eulers[2]) )
    reconstructedODF[bin] += unitInc

  return orientations, reconstructedODF


def TothVanHoutteSTAT (ODF,nSamples):
  """ ODF contains 'dV_V' (normalized to 1), 'center', 'intervals', 'limits' (in radians) """

  nBins = ODF['intervals'][0]*ODF['intervals'][1]*ODF['intervals'][2]
  deltas = [limit/intervals for limit,intervals in zip(ODF['limits'],ODF['intervals'])]
  orientations = [None]*nSamples
  reconstructedODF = [0.0]*nBins
  unitInc = 1.0/nSamples
  
  selectors = [random.random() for i in range(nSamples)]
  selectors.sort()
  indexSelector = 0
  
  cumdV_V = 0.0
  countSamples = 0
  
  for bin in range(nBins) :
    cumdV_V += ODF['dV_V'][bin]
    while indexSelector < nSamples and selectors[indexSelector] < cumdV_V:
      Eulers = binAsEulers(bin,ODF['intervals'],deltas,ODF['center'])
      orientations[countSamples] = '%g\t%g\t%g' %( math.degrees(Eulers[0]),math.degrees(Eulers[1]),math.degrees(Eulers[2]) )
      reconstructedODF[bin] += unitInc
      countSamples += 1
      indexSelector += 1
      
  print 'created set of',countSamples,'when asked to deliver',nSamples
  
  return orientations, reconstructedODF


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Transform linear binned data into Euler angles.
""", version=string.replace(scriptID,'\n','\\n')
)

parser.add_option('-f', '--file', dest='file', type='string', metavar = 'string', \
                     help='file name, the input file is generated by the script "OIMlinear2linearODF.py"')
parser.add_option('-o', '--output', dest='output', type='string', metavar = 'string', \
                     help='the prefix of output files name.')
parser.add_option('-n', '--number', dest='number', type='int', metavar = 'int', \
                     help='the number of orientations needed to be generated, the default is [%default]')
parser.add_option('-a','--algorithm', dest='algorithm', type='string', metavar = 'string', \
                     help='''The algorithm adopted, three algorithms are provided, 
                           that is:                      
                             [IA]:   direct inversion,                        
                             [STAT]: Van Houtte,                 
                             [MC]:   Monte Carlo.                      
                           the default is [%default].''')
parser.add_option('-p','--phase', dest='phase', type='int', metavar = 'int',
                  help='phase index to be used [%default]')
parser.add_option('--crystallite', dest='crystallite', type='int', metavar = 'int',
                  help='crystallite index to be used [%default]')
parser.add_option('-c', '--configuration', dest='config', action='store_true',
                  help='output material configuration [%default]')

parser.set_defaults(number      = 500)
parser.set_defaults(output      = 'texture')
parser.set_defaults(algorithm   = 'IA')
parser.set_defaults(phase       = 1)
parser.set_defaults(crystallite = 1)

options = parser.parse_args()[0]
# check usage
if not os.path.exists(options.file):
  parser.error('binnedODF file does not exist'); sys.exit()

nameBinnedODF  = options.file
nameSampledODF = options.output
nSamples       = options.number
methods        = [options.algorithm]

#open binned ODF
try:
  fileBinnedODF = open(nameBinnedODF,'r')
except:
  print 'unable to open binnedODF:', nameBinnedODF;
  sys.exit(1);

# process header info
ODF = {}
ODF['limits'] = [math.radians(float(limit)) for limit in fileBinnedODF.readline().split()[0:3]]
ODF['deltas'] = [math.radians(float(delta)) for delta in fileBinnedODF.readline().split()[0:3]]
ODF['intervals'] = [int(interval) for interval in fileBinnedODF.readline().split()[0:3]]
#ODF['intervals'] = [int(round(limit/delta)) for limit,delta in zip(ODF['limits'],ODF['deltas'])]
nBins = ODF['intervals'][0]*ODF['intervals'][1]*ODF['intervals'][2]

print 'Limit:    ', [math.degrees(limit) for limit in ODF['limits']]
print 'Delta:    ', [math.degrees(delta) for delta in ODF['deltas']]
print 'Interval: ', ODF['intervals']

centering = fileBinnedODF.readline()
if re.search('cell',centering.lower()):
  ODF['center'] = 0.5
  print 'cell-centered data (offset %g)'%ODF['center']
else:
  ODF['center'] = 0.0
  print 'vertex-centered data (offset %g)'%ODF['center']

fileBinnedODF.readline() # skip blank delimiter

# read linear binned data
linesBinnedODF = fileBinnedODF.readlines()
fileBinnedODF.close()

if len(linesBinnedODF) != nBins:
  print 'expecting', nBins, 'values but got', len(linesBinnedODF)
  sys.exit(1)

# build binnedODF array
print 'reading',nBins,'values'
sumdV_V = 0.0
ODF['dV_V'] = [None]*nBins
ODF['nNonZero'] = 0
dg = ODF['deltas'][0]*2*math.sin(ODF['deltas'][1]/2.0)*ODF['deltas'][2]
for bin in range(nBins):
  ODF['dV_V'][bin] = \
  max(0.0,float(linesBinnedODF[bin])) * dg * \
  math.sin(((bin//ODF['intervals'][2])%ODF['intervals'][1]+ODF['center'])*ODF['deltas'][1])
  if ODF['dV_V'][bin] > 0.0:
    sumdV_V += ODF['dV_V'][bin]
    ODF['nNonZero'] += 1

for bin in range(nBins): ODF['dV_V'][bin] /= sumdV_V # normalize dV/V

print 'non-zero fraction:', float(ODF['nNonZero'])/nBins,'(%i/%i)'%(ODF['nNonZero'],nBins)
print 'Volume integral of ODF:', sumdV_V
print 'Reference Integral:', ODF['limits'][0]*ODF['limits'][2]*(1-math.cos(ODF['limits'][1]))

# call methods
Functions = {'IA': 'directInversion', 'STAT': 'TothVanHoutteSTAT', 'MC': 'MonteCarloBins'}
Orientations = {}
ReconstructedODF = {}
for method in methods:
  Orientations[method], ReconstructedODF[method] = (globals()[Functions[method]])(ODF,nSamples)

# calculate accuracy of sample
squaredDiff = {}
squaredRelDiff = {}
mutualProd = {}
indivSum = {}
indivSquaredSum = {}
for method in ['orig']+methods:
  squaredDiff[method] = 0.0
  squaredRelDiff[method] = 0.0
  mutualProd[method] = 0.0
  indivSum[method] = 0.0
  indivSquaredSum[method] = 0.0

for bin in range(nBins):
  for method in methods:
    squaredDiff[method] += (ODF['dV_V'][bin] - ReconstructedODF[method][bin])**2
    if ODF['dV_V'][bin] > 0.0:
      squaredRelDiff[method] += (ODF['dV_V'][bin] - ReconstructedODF[method][bin])**2/ODF['dV_V'][bin]**2
    mutualProd[method] += ODF['dV_V'][bin]*ReconstructedODF[method][bin]
    indivSum[method] += ReconstructedODF[method][bin]
    indivSquaredSum[method] += ReconstructedODF[method][bin]**2
  indivSum['orig'] += ODF['dV_V'][bin]
  indivSquaredSum['orig'] += ODF['dV_V'][bin]**2

print 'sqrt(N*)RMSD of ODFs:\t', [math.sqrt(nSamples*squaredDiff[method]) for method in methods]
print 'RMSrD of ODFs:\t', [math.sqrt(squaredRelDiff[method]) for method in methods]
print 'rMSD of ODFs:\t', [squaredDiff[method]/indivSquaredSum['orig'] for method in methods]
#print 'correlation slope:\t', [(nBins*mutualProd[method]-indivSum['orig']*indivSum[method])/(nBins*indivSquaredSum['orig']-indivSum['orig']**2) for method in ('IA','STAT','MC')]
#print 'correlation confidence:\t', [(mutualProd[method]-indivSum['orig']*indivSum[method]/nBins)/\
#                  (nBins*math.sqrt((indivSquaredSum['orig']/nBins-(indivSum['orig']/nBins)**2)*(indivSquaredSum[method]/nBins-(indivSum[method]/nBins)**2))) for method in ('IA','STAT','MC')]
print 'nNonZero correlation slope:\t', [(ODF['nNonZero']*mutualProd[method]-indivSum['orig']*indivSum[method])/(ODF['nNonZero']*indivSquaredSum['orig']-indivSum['orig']**2) for method in methods]
print 'nNonZero correlation confidence:\t', [(mutualProd[method]-indivSum['orig']*indivSum[method]/ODF['nNonZero'])/\
                  (ODF['nNonZero']*math.sqrt((indivSquaredSum['orig']/ODF['nNonZero']-(indivSum['orig']/ODF['nNonZero'])**2)*(indivSquaredSum[method]/ODF['nNonZero']-(indivSum[method]/ODF['nNonZero'])**2))) for method in methods]

for method in  methods:
  if method == 'IA' and nSamples < ODF['nNonZero']:
    strOpt = '(%i)'%ODF['nNonZero']
  else:
    strOpt = ''
  try:
    fileSampledODF = open(nameSampledODF+'.'+method+'sampled_'+str(nSamples)+strOpt, 'w')
    fileSampledODF.write('%i\n'%nSamples)
    fileSampledODF.write('\n'.join(Orientations[method])+'\n')
    fileSampledODF.close()
  except:
    print 'unable to write sampledODF:', nameSampledODF+'.'+method+'sampled_'+str(nSamples)+strOpt

  try:
    fileRegressionODF = open(nameSampledODF+'.'+method+'regression_'+str(nSamples)+strOpt, 'w')
    fileRegressionODF.write('\n'.join([a+'\t'+b for (a,b) in zip(map(str,ReconstructedODF[method]),map(str,ODF['dV_V']))])+'\n')
    fileRegressionODF.close()
  except:
    print 'unable to write RegressionODF:', nameSampledODF+'.'+method+'regression_'+str(nSamples)+strOpt
  
  if options.config:
    try: 
      fileConfig = open(nameSampledODF+'.'+method+str(nSamples)+'.config', 'w')                                                                           # write config file
      formatwidth = 1
      fileConfig.write('#-------------------#')
      fileConfig.write('\n<microstructure>\n')
      fileConfig.write('#-------------------#\n')

      for i,ID in enumerate(xrange(nSamples)):
        fileConfig.write('[Grain%s]\n'%(str(ID+1).zfill(formatwidth)) + \
                         'crystallite %i\n'%options.crystallite + \
                         '(constituent)   phase %i   texture %s   fraction 1.0\n'%(options.phase,str(ID+1).rjust(formatwidth)))

      fileConfig.write('\n#-------------------#')
      fileConfig.write('\n<texture>\n')
      fileConfig.write('#-------------------#\n')
      for ID in xrange(nSamples):
        eulers = re.split(r'[\t]', Orientations[method][ID].strip())

        fileConfig.write('[Grain%s]\n'%(str(ID+1).zfill(formatwidth)) + \
                         '(gauss)   phi1 %10.5f   Phi %10.5f   phi2 %10.6f   scatter 0.0   fraction 1.0\n'\
                         %(float(eulers[0]),float(eulers[1]),float(eulers[2])))
      fileConfig.close()
    except:
      print 'unable to write material.config file:', nameSampledODF+'.'+method+str(nSamples)+'.config'
