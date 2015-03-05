#!/usr/bin/env python

import os,sys,re

try:
  inName = sys.argv[1]
  outName = os.path.splitext(inName)[0]+'.linearODF'
  nPhi1,nPHI,nPhi2 = map(int,sys.argv[2:5])
except:
  print "\nusage:",sys.argv[0],"file nPhi1 nPHI nPhi2\n"
  sys.exit(1)

N = (nPhi1-1)*(nPHI-1)*(nPhi2-1)

try:
  inFile = open(inName,'r')
  content = inFile.readlines()
except:
  print 'unable to read:',inName
  sys.exit(1)
try:
  outFile = open(outName,'w')
except:
  print 'unable to write:',outName
  sys.exit(1)

ODF = [[[[None] for k in range(nPhi2)] for j in range(nPHI)] for i in range(nPhi1)]
linear = [None]*N
line = 0

while (content[line].startswith('#')):  # skip comments at start of file
  line += 1

for iPhi1 in range(nPhi1):
  for iPHI in range(nPHI):
    for iPhi2 in range(nPhi2):
      words = content[line].split()
      ODF[iPhi1][iPHI][iPhi2] = float(words[3]) # extract intensity (in column 4)
      line += 1

for iPhi1 in range(nPhi1-1):
  for iPHI in range(nPHI-1):
    for iPhi2 in range(nPhi2-1):
      linear[iPhi1*(nPHI-1)*(nPhi2-1)+iPHI*(nPhi2-1)+iPhi2] = (\
        ODF[iPhi1  ][iPHI  ][iPhi2  ] + \
        ODF[iPhi1  ][iPHI  ][iPhi2+1] + \
        ODF[iPhi1  ][iPHI+1][iPhi2  ] + \
        ODF[iPhi1  ][iPHI+1][iPhi2+1] + \
        ODF[iPhi1+1][iPHI  ][iPhi2  ] + \
        ODF[iPhi1+1][iPHI  ][iPhi2+1] + \
        ODF[iPhi1+1][iPHI+1][iPhi2  ] + \
        ODF[iPhi1+1][iPHI+1][iPhi2+1] \
      ) / 8.0


inFile.close()

outFile.write('rangePhi1?\trangePHI?\trangePhi2?\n')
outFile.write('%i\t%i\t%i needs to be converted to angular steps\n'%(nPhi1-1,nPHI-1,nPhi2-1))
outFile.write('cell-centered data\n')
outFile.write('\n')

for i in range(N):
  outFile.write('%g\n'%(linear[i]))
  
outFile.close()

