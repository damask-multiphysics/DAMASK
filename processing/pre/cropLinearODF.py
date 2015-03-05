#!/usr/bin/env python

import os,sys,math,re

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



# check usage
try:
  inName = sys.argv[1]
  outLimits = sys.argv[2:5]
except:
  print "usage:",sys.argv[0],"nameLinearODF limitPhi1 limitPHI limitPhi2"
  sys.exit(1);


#open binned ODF
try:
  inFile = open(inName,'r')
except:
  print 'unable to open binnedODF:', inName;
  sys.exit(1);

# process header info
ODF = {}
inLimits = [math.radians(int(float(limit))) for limit in inFile.readline().split()]
outLimits = [math.radians(int(float(limit))) for limit in outLimits]
deltas = [math.radians(float(delta)) for delta in inFile.readline().split()]
inIntervals = [int(limit/delta) for limit,delta in zip(inLimits,deltas)]
outIntervals = [int(limit/delta) for limit,delta in zip(outLimits,deltas)]
inBins = inIntervals[0]*inIntervals[1]*inIntervals[2]

print 'Limit:', [math.degrees(limit) for limit in inLimits]
print 'Crop:', [math.degrees(limit) for limit in outLimits]
print 'Delta:', [math.degrees(delta) for delta in deltas]
print 'Interval:', inIntervals
print 'Interval:', outIntervals

centering = inFile.readline()
if re.search('cell',centering.lower()):
  ODF['center'] = 0.5
  print 'cell-centered data (offset %g)'%ODF['center']
else:
  ODF['center'] = 0.0
  print 'vertex-centered data (offset %g)'%ODF['center']

inFile.readline() # skip blank delimiter

# read linear binned data
inODF = map(float,inFile.readlines())
inFile.close()

if len(inODF) != inBins:
  print 'expecting', inBins, 'values but got', len(inODF)
  sys.exit(1)

try:
  outName = os.path.splitext(inName)[0]+'_%ix%ix%i'%(outIntervals[0],outIntervals[1],outIntervals[2])+'.linearODF'
  outFile = open(outName,'w')
except:
  print 'unable to write:',outName
  sys.exit(1)

outFile.write('%g\t%g\t%g\n'%(\
    math.degrees(outIntervals[0]*deltas[0]),\
    math.degrees(outIntervals[1]*deltas[1]),\
    math.degrees(outIntervals[2]*deltas[2]) ))
outFile.write('%i\t%i\t%i\n'%(math.degrees(deltas[0]),math.degrees(deltas[1]),math.degrees(deltas[2])))
outFile.write('%s-centered data\n'%{True:'vertex',False:'cell'}[ODF['center']==0.0])
outFile.write('\n')

for phi1 in range(outIntervals[0]):
  for Phi in range(outIntervals[1]):
    for phi2 in range(outIntervals[2]):
      outFile.write('%g\n'%(inODF[((phi1*inIntervals[1])+Phi)*inIntervals[2]+phi2]))

outFile.close()
