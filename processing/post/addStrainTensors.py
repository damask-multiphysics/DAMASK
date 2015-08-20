#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

def operator(stretch,strain,eigenvalues):
  ''' Albrecht Bertram: Elasticity and Plasticity of Large Deformations An Introduction (3rd Edition, 2012), p. 102 '''
  return {
    'V#ln':    np.log(eigenvalues)                                 ,
    'U#ln':    np.log(eigenvalues)                                 ,
    'V#Biot':  ( np.ones(3,'d') - 1.0/eigenvalues )                ,
    'U#Biot':  ( eigenvalues - np.ones(3,'d') )                    ,
    'V#Green': ( np.ones(3,'d') - 1.0/eigenvalues/eigenvalues) *0.5,
    'U#Green': ( eigenvalues*eigenvalues - np.ones(3,'d'))     *0.5,
         }[stretch+'#'+strain]


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing given strains based on given stretches of requested deformation gradient column(s).

""", version = scriptID)

parser.add_option('-u','--right',
                  dest = 'right',
                  action = 'store_true',
                  help = 'material strains based on right Cauchy--Green deformation, i.e., C and U')
parser.add_option('-v','--left',
                  dest = 'left',
                  action = 'store_true',
                  help = 'spatial strains based on left Cauchy--Green deformation, i.e., B and V')
parser.add_option('-0','--logarithmic',
                  dest = 'logarithmic',
                  action = 'store_true',
                  help = 'calculate logarithmic strain tensor')
parser.add_option('-1','--biot',
                  dest = 'biot',
                  action = 'store_true',
                  help = 'calculate biot strain tensor')
parser.add_option('-2','--green',
                  dest = 'green',
                  action = 'store_true',
                  help = 'calculate green strain tensor')
parser.add_option('-f','--defgrad',
                  dest = 'defgrad',
                  action = 'extend',
                  metavar = '<string LIST>',
                  help = 'heading(s) of columns containing deformation tensor values [%default]')

parser.set_defaults(right       = False,
                    left        = False,
                    logarithmic = False,
                    biot        = False,
                    green       = False,
                    defgrad     = ['f'],
                   )

(options,filenames) = parser.parse_args()

stretches = []
strains = []

if options.right: stretches.append('U')
if options.left:  stretches.append('V')
if options.logarithmic: strains.append('ln')
if options.biot:        strains.append('Biot')
if options.green:       strains.append('Green')

if options.defgrad == None:
  parser.error('no data column specified.')

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except: continue
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name else ''))

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  items = {
            'tensor': {'dim': 9, 'shape': [3,3], 'labels':options.defgrad, 'column': []},
          }
  errors  = []
  remarks = []
  
  for type, data in items.iteritems():
    for what in data['labels']:
      dim = table.label_dimension(what)
      if dim != data['dim']: remarks.append('column {} is not a {}...'.format(what,type))
      else:
        items[type]['column'].append(table.label_index(what))
        for theStretch in stretches:
          for theStrain in strains:
            table.labels_append(['{}_{}({}){}'.format(i+1,                                          # extend ASCII header with new labels  
                                                      theStrain,
                                                      theStretch,
                                                      what if what != 'f' else '') for i in xrange(9)])

  if remarks != []: table.croak(remarks)
  if errors  != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  stretch = {}
  outputAlive = True

  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for column in items['tensor']['column']:                                                        # loop over all requested defgrads
      F = np.array(map(float,table.data[column:column+items['tensor']['dim']]),'d').reshape(items['tensor']['shape'])
      (U,S,Vh) = np.linalg.svd(F)                                                                   # singular value decomposition
      R = np.dot(U,Vh)                                                                              # rotation of polar decomposition
      stretch['U'] = np.dot(np.linalg.inv(R),F)                                                     # F = RU
      stretch['V'] = np.dot(F,np.linalg.inv(R))                                                     # F = VR

      for theStretch in stretches:
        stretch[theStretch] = np.where(abs(stretch[theStretch]) < 1e-12, 0, stretch[theStretch])    # kill nasty noisy data
        (D,V) = np.linalg.eig(stretch[theStretch])                                                  # eigen decomposition (of symmetric matrix)
        neg = np.where(D < 0.0)                                                                     # find negative eigenvalues ...
        D[neg]   *= -1.                                                                             # ... flip value ...
        V[:,neg] *= -1.                                                                             # ... and vector
        for i,eigval in enumerate(D):
          if np.dot(V[:,i],V[:,(i+1)%3]) != 0.0:                                                    # check each vector for orthogonality
              V[:,(i+1)%3] = np.cross(V[:,(i+2)%3],V[:,i])                                          # correct next vector
              V[:,(i+1)%3] /= np.sqrt(np.dot(V[:,(i+1)%3],V[:,(i+1)%3].conj()))                     # and renormalize (hyperphobic?)
        for theStrain in strains:
          d = operator(theStretch,theStrain,D)                                                      # operate on eigenvalues of U or V
          eps = (np.dot(V,np.dot(np.diag(d),V.T)).real).reshape(9)                                  # build tensor back from eigenvalue/vector basis

          table.data_append(list(eps))

# ------------------------------------------ output result -----------------------------------------

    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
