#!/usr/bin/env python

from optparse import OptionParser
import damask
import os,sys,math,re,random,string
import numpy as np

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = scriptID.split()[1]

random.seed(1)

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
  return [int((euler+(0.5-center)*delta)//delta)%interval \
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


# ---------------------- sampling methods -----------------------------------------------------------------------

# ----- efficient algorithm ---------

def directInversion (ODF,nSamples):
  """ ODF contains 'dV_V' (normalized to 1), 'center', 'intervals', 'limits' (in radians) """
  
  nOptSamples = max(ODF['nNonZero'],nSamples)         # random subsampling if too little samples requested

  nInvSamples = 0
  repetition = [None]*ODF['nBins']
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
  repetition = [None]*ODF['nBins'] # preallocate and clear
    
  for bin in range(ODF['nBins']): # loop over bins
    repetition[bin] = int(round(ODF['dV_V'][bin]*scale)) # calc repetition

  # build set
  set = [None]*nInvSamples
  i = 0
  for bin in range(ODF['nBins']):
    set[i:i+repetition[bin]] = [bin]*repetition[bin] # fill set with bin, i.e. orientation
    i += repetition[bin] # advance set counter
  
  orientations = [None]*nSamples
  reconstructedODF = [0.0]*ODF['nBins']
  unitInc = 1.0/nSamples
  for j in range(nSamples):
    if (j == nInvSamples-1): ex = j
    else: ex = int(round(random.uniform(j+0.5,nInvSamples-0.5)))
    bin = set[ex]
    bins = binAsBins(bin,ODF['interval'])
    Eulers = binAsEulers(bin,ODF['interval'],ODF['delta'],ODF['center'])
    orientations[j] = '%g\t%g\t%g' %( math.degrees(Eulers[0]),math.degrees(Eulers[1]),math.degrees(Eulers[2]) )
    reconstructedODF[bin] += unitInc
    set[ex] = set[j] # exchange orientations
  
  return orientations, reconstructedODF


# ----- trial and error algorithms ---------

def MonteCarloEulers (ODF,nSamples):
  """ ODF contains 'dV_V' (normalized to 1), 'center', 'intervals', 'limits' (in radians) """
  
  countMC = 0
  maxdV_V = max(ODF['dV_V'])
  orientations = [None]*nSamples
  reconstructedODF = [0.0]*ODF['nBins']
  unitInc = 1.0/nSamples
  
  for j in range(nSamples):
    MC = maxdV_V*2.0
    bin = 0
    while MC > ODF['dV_V'][bin]:
      countMC += 1
      MC = maxdV_V*random.random()
    Eulers = [limit*random.random() for limit in ODF['limit']]
    bins = EulersAsBins(Eulers,ODF['interval'],ODF['delta'],ODF['center'])
    bin = binsAsBin(bins,ODF['interval'])
    orientations[j] = '%g\t%g\t%g' %( math.degrees(Eulers[0]),math.degrees(Eulers[1]),math.degrees(Eulers[2]) )
    reconstructedODF[bin] += unitInc

  return orientations, reconstructedODF, countMC


def MonteCarloBins (ODF,nSamples):
  """ ODF contains 'dV_V' (normalized to 1), 'center', 'intervals', 'limits' (in radians) """
  
  countMC = 0
  maxdV_V = max(ODF['dV_V'])
  orientations = [None]*nSamples
  reconstructedODF = [0.0]*ODF['nBins']
  unitInc = 1.0/nSamples
  
  for j in range(nSamples):
    MC = maxdV_V*2.0
    bin = 0
    while MC > ODF['dV_V'][bin]:
      countMC += 1
      MC  = maxdV_V*random.random()
      bin = int(ODF['nBins'] * random.random())
    Eulers = binAsEulers(bin,ODF['interval'],ODF['delta'],ODF['center'])
    orientations[j] = '%g\t%g\t%g' %( math.degrees(Eulers[0]),math.degrees(Eulers[1]),math.degrees(Eulers[2]) )
    reconstructedODF[bin] += unitInc

  return orientations, reconstructedODF


def TothVanHoutteSTAT (ODF,nSamples):
  """ ODF contains 'dV_V' (normalized to 1), 'center', 'intervals', 'limits' (in radians) """

  orientations = [None]*nSamples
  reconstructedODF = [0.0]*ODF['nBins']
  unitInc = 1.0/nSamples
  
  selectors = [random.random() for i in range(nSamples)]
  selectors.sort()
  indexSelector = 0
  
  cumdV_V = 0.0
  countSamples = 0
  
  for bin in range(ODF['nBins']) :
    cumdV_V += ODF['dV_V'][bin]
    while indexSelector < nSamples and selectors[indexSelector] < cumdV_V:
      Eulers = binAsEulers(bin,ODF['interval'],ODF['delta'],ODF['center'])
      orientations[countSamples] = '%g\t%g\t%g' %( math.degrees(Eulers[0]),math.degrees(Eulers[1]),math.degrees(Eulers[2]) )
      reconstructedODF[bin] += unitInc
      countSamples += 1
      indexSelector += 1
      
  print 'created set of',countSamples,'when asked to deliver',nSamples
  
  return orientations, reconstructedODF


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
identifiers = {
        'limit':   ['phi1','phi','phi2'],
        'delta':   ['phi1','phi','phi2'],
          }
mappings = {
        'limit':       lambda x: math.radians(float(x)),
        'delta':       lambda x: math.radians(float(x)),
        'origin':      lambda x: str(x),
          }

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Transform linear binned data into Euler angles.

""", version = scriptID)

parser.add_option('-n', '--number', dest='number', type='int', metavar = 'int',
                     help='the number of orientations needed to be generated, the default is [%default]')
parser.add_option('-a','--algorithm', dest='algorithm', type='string', metavar = 'string',
                     help='''The algorithm adopted, three algorithms are provided, 
                           that is:                      
                             [IA]:   direct inversion,                        
                             [STAT]: Van Houtte,                 
                             [MC]:   Monte Carlo.                      
                           the default is [%default].''') #make (multiple) choice
parser.add_option('-p','--phase', dest='phase', type='int', metavar = 'int',
                  help='phase index to be used [%default]')
parser.add_option('--crystallite', dest='crystallite', type='int', metavar = 'int',
                  help='crystallite index to be used [%default]')

parser.set_defaults(number      = 500)
parser.set_defaults(algorithm   = 'IA')
parser.set_defaults(phase       = 1)
parser.set_defaults(crystallite = 1)
(options,filenames) = parser.parse_args()

nSamples       = options.number
methods        = [options.algorithm]

#--- setup file handles ---------------------------------------------------------------------------
files = []
if filenames == []:
  files.append({'name':'STDIN',
                'input':sys.stdin,
                'output':sys.stdout,
                'croak':sys.stderr,
               })
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name,
                    'input':open(name),
                    'output':open(name+'_tmp','w'),
                    'croak':sys.stdout,
                    })

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  file['croak'].write('\033[1m' + scriptName + '\033[0m: ' + (file['name'] if file['name'] != 'STDIN' else '') + '\n')

  ODF = {
          'limit':   np.empty(3,'d'),
          'delta':   np.empty(3,'d'),
          'interval':np.empty(3,'i'),
          'origin':  '' 
         }

  table = damask.ASCIItable(file['input'],file['output'],buffered = False)
  table.head_read()

  for header in table.info:
    headitems = map(str.lower,header.split())
    if len(headitems) == 0: continue
    if headitems[0] in mappings.keys():
      if headitems[0] in identifiers.keys():
        for i in xrange(len(identifiers[headitems[0]])):
          ODF[headitems[0]][i] = \
            mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])
      else:
        ODF[headitems[0]] = mappings[headitems[0]](headitems[1])

  ODF['interval'] = np.array([int(round(limit/delta)) for limit,delta in zip(ODF['limit'],ODF['delta'])])
  ODF['nBins'] = ODF['interval'].prod()
  
  if re.search('boundary',ODF['origin'].lower()):
    ODF['center'] = 0.5
  else:
    ODF['center'] = 0.0
  
  binnedODF=table.data_readArray([table.labels.index('intensity')])

  if binnedODF[0] != ODF['nBins']:
    print 'expecting', ODF['nBins'], 'values but got', len(linesBinnedODF)
    sys.exit(1)
  
  # build binnedODF array
  sumdV_V = 0.0
  ODF['dV_V'] = [None]*ODF['nBins']
  ODF['nNonZero'] = 0
  dg = ODF['delta'][0]*2*math.sin(ODF['delta'][1]/2.0)*ODF['delta'][2]
  for bin in range(ODF['nBins']):
    ODF['dV_V'][bin] = \
    max(0.0,table.data[bin,0]) * dg * \
    math.sin(((bin//ODF['interval'][2])%ODF['interval'][1]+ODF['center'])*ODF['delta'][1])
    if ODF['dV_V'][bin] > 0.0:
      sumdV_V += ODF['dV_V'][bin]
      ODF['nNonZero'] += 1
  
  for bin in range(ODF['nBins']): ODF['dV_V'][bin] /= sumdV_V # normalize dV/V
  
  print 'non-zero fraction:', float(ODF['nNonZero'])/ODF['nBins'],'(%i/%i)'%(ODF['nNonZero'],ODF['nBins'])
  print 'Volume integral of ODF:', sumdV_V
  print 'Reference Integral:', ODF['limit'][0]*ODF['limit'][2]*(1-math.cos(ODF['limit'][1]))
  
  # call methods
  Functions = {'IA': 'directInversion', 'STAT': 'TothVanHoutteSTAT', 'MC': 'MonteCarloBins'}
  method = Functions[options.algorithm]

  Orientations, ReconstructedODF = (globals()[method])(ODF,nSamples)
  
  # calculate accuracy of sample
  squaredDiff     = {'orig':0.0,method:0.0}
  squaredRelDiff  = {'orig':0.0,method:0.0}
  mutualProd      = {'orig':0.0,method:0.0}
  indivSum        = {'orig':0.0,method:0.0}
  indivSquaredSum = {'orig':0.0,method:0.0}
  
  for bin in range(ODF['nBins']):
    squaredDiff[method] += (ODF['dV_V'][bin] - ReconstructedODF[bin])**2
    if ODF['dV_V'][bin] > 0.0:
      squaredRelDiff[method] += (ODF['dV_V'][bin] - ReconstructedODF[bin])**2/ODF['dV_V'][bin]**2
    mutualProd[method] += ODF['dV_V'][bin]*ReconstructedODF[bin]
    indivSum[method] += ReconstructedODF[bin]
    indivSquaredSum[method] += ReconstructedODF[bin]**2
  indivSum['orig'] += ODF['dV_V'][bin]
  indivSquaredSum['orig'] += ODF['dV_V'][bin]**2
  
  print 'sqrt(N*)RMSD of ODFs:\t', math.sqrt(nSamples*squaredDiff[method])
  print 'RMSrD of ODFs:\t',        math.sqrt(squaredRelDiff[method])
  print 'rMSD of ODFs:\t',         squaredDiff[method]/indivSquaredSum['orig']
  print 'nNonZero correlation slope:\t', (ODF['nNonZero']*mutualProd[method]-indivSum['orig']*indivSum[method])/\
                                         (ODF['nNonZero']*indivSquaredSum['orig']-indivSum['orig']**2)
  print 'nNonZero correlation confidence:\t',\
            (mutualProd[method]-indivSum['orig']*indivSum[method]/ODF['nNonZero'])/\
             (ODF['nNonZero']*math.sqrt((indivSquaredSum['orig']/ODF['nNonZero']-(indivSum['orig']/ODF['nNonZero'])**2)*\
                                        (indivSquaredSum[method]/ODF['nNonZero']-(indivSum[method]/ODF['nNonZero'])**2)))
  
  if method == 'IA' and nSamples < ODF['nNonZero']:
    strOpt = '(%i)'%ODF['nNonZero']
  
  formatwidth = 1
  file['output'].write('#-------------------#')
  file['output'].write('\n<microstructure>\n')
  file['output'].write('#-------------------#\n')
  
  for i,ID in enumerate(xrange(nSamples)):
    file['output'].write('[Grain%s]\n'%(str(ID+1).zfill(formatwidth)) + \
                     'crystallite %i\n'%options.crystallite + \
                     '(constituent)   phase %i   texture %s   fraction 1.0\n'%(options.phase,str(ID+1).rjust(formatwidth)))
  
  file['output'].write('\n#-------------------#')
  file['output'].write('\n<texture>\n')
  file['output'].write('#-------------------#\n')
  for ID in xrange(nSamples):
    eulers = re.split(r'[\t]', Orientations[ID].strip())
  
    file['output'].write('[Grain%s]\n'%(str(ID+1).zfill(formatwidth)) + \
                     '(gauss)   phi1 %10.5f   Phi %10.5f   phi2 %10.6f   scatter 0.0   fraction 1.0\n'\
                     %(float(eulers[0]),float(eulers[1]),float(eulers[2])))
  #--- output finalization -------------------------------------------------------------------------- 
  if file['name'] != 'STDIN':
     file['output'].close()
     os.rename(file['name']+'_tmp',
            os.path.splitext(file['name'])[0] +'_'+method+'%s'%('_material.config'))
