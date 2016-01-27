#!/usr/bin/python
# -*- coding: UTF-8 no BOM -*-

import os,string,sys,re
from optparse import OptionParser
import numpy as np
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

sampleSym = {  'Orthotropic'  : (90,90,90),
               'Triclinic'    : (360,180,360)
            }

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Transform the binned texture data from "TSL OIM Analysis" into linear ODF data,

""", version = scriptID)

parser.add_option('-s', '--symmetry', dest='symmetry', choices=sampleSym.keys(),
                     help='Sample symmetry {%s} [Triclinic]'%(' '.join(sampleSym.keys())))
                                     
parser.set_defaults(symmetry = 'Triclinic')

(options,filenames) = parser.parse_args()

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

  while True:                                                                                       # read header (forward and get bin Size)
    line = file['input'].readline()
    words = line.split()
    if len(words)>=3:
      if words[1]=='Bin' and words[2]=='Size:': binSize=float(words[3][:-1])
    if not line.startswith('#'): break

  delta = [sampleSym[options.symmetry][i]/binSize for i in xrange(3)]

  nPhi1,nPHI,nPhi2 = map(int,delta)
  dPhi1,dPHI,dPhi2 = [sampleSym[options.symmetry][i]/delta[i] for i in xrange(3)]

  N = (nPhi1-1)*(nPHI-1)*(nPhi2-1)


  ODF = [[[[None] for k in range(nPhi2)] for j in range(nPHI)] for i in range(nPhi1)]
  linear = [None]*N

  ODF = np.empty([nPhi1,nPHI,nPhi2],'d')

  for iPhi1 in range(nPhi1):
    for iPHI in range(nPHI):
      for iPhi2 in range(nPhi2):
        ODF[iPhi1,iPHI,iPhi2] = float(line.split()[3])*0.125                                      # extract intensity (in column 4) and weight by 1/8 (since we convert from the 8 corners to the center later on)
        line = file['input'].readline()

  for iPhi1 in range(nPhi1-1):
    for iPHI in range(nPHI-1):
      for iPhi2 in range(nPhi2-1):
        linear[iPhi1*(nPHI-1)*(nPhi2-1)+iPHI*(nPhi2-1)+iPhi2] =\
          ODF[iPhi1  ,iPHI  ,iPhi2  ] +\
          ODF[iPhi1  ,iPHI  ,iPhi2+1] +\
          ODF[iPhi1  ,iPHI+1,iPhi2  ] +\
          ODF[iPhi1  ,iPHI+1,iPhi2+1] +\
          ODF[iPhi1+1,iPHI  ,iPhi2  ] +\
          ODF[iPhi1+1,iPHI  ,iPhi2+1] +\
          ODF[iPhi1+1,iPHI+1,iPhi2  ] +\
          ODF[iPhi1+1,iPHI+1,iPhi2+1]


  file['output'].write('4 header\n')
  file['output'].write('limit phi1 %-6.2f Phi %-6.2f phi2 %-6.2f\n'%sampleSym[options.symmetry])
  file['output'].write('delta phi1 %-6.2f Phi %-6.2f phi2 %-6.2f\n'%(dPhi1,dPHI,dPhi2))
  file['output'].write('centration cell-centered\n')
  file['output'].write('density\n')

  for i in range(N):
    file['output'].write('%g\n'%(linear[i]))
  
#--- output finalization -------------------------------------------------------------------------- 
  if file['name'] != 'STDIN':
     file['output'].close()
     os.rename(file['name']+'_tmp',os.path.splitext(file['name'])[0] +'.linearODF')
