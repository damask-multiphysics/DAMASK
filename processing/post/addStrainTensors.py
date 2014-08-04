#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,string
import numpy as np
from collections import defaultdict
from optparse import OptionParser
import damask

scriptID = '$Id$'
scriptName = scriptID.split()[1]

def operator(stretch,strain,eigenvalues):
  return {
    'V#ln':    np.log(eigenvalues)                                 ,
    'U#ln':    np.log(eigenvalues)                                 ,
    'V#Biot':  ( np.ones(3,'d') - 1.0/eigenvalues )                ,
    'U#Biot':  ( eigenvalues - np.ones(3,'d') )                    ,
    'V#Green': ( np.ones(3,'d') - 1.0/eigenvalues*eigenvalues) *0.5,
    'U#Green': ( eigenvalues*eigenvalues - np.ones(3,'d'))     *0.5,
         }[stretch+'#'+strain]


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing given strains based on given stretches of requested deformation gradient column(s).

""" + string.replace(scriptID,'\n','\\n')
)

parser.add_option('-u','--right',       dest='right', action='store_true',
                                        help='material strains based on right Cauchy--Green deformation, i.e., C and U [%default]')
parser.add_option('-v','--left',        dest='left', action='store_true',
                                        help='spatial strains based on left Cauchy--Green deformation, i.e., B and V [%default]')
parser.add_option('-l','-0','--logarithmic', dest='logarithmic', action='store_true',
                                        help='calculate logarithmic strain tensor [%default]')
parser.add_option('-b','-1','--biot',   dest='biot', action='store_true',
                                        help='calculate biot strain tensor [%default]')
parser.add_option('-g','-2','--green',  dest='green', action='store_true',
                                        help='calculate green strain tensor [%default]')
parser.add_option('-f','--defgrad',     dest='defgrad', action='extend', type='string', metavar = '<string LIST>',
                                        help='heading(s) of columns containing deformation tensor values %default')
parser.set_defaults(right       = False)
parser.set_defaults(left        = False)
parser.set_defaults(logarithmic = False)
parser.set_defaults(biot        = False)
parser.set_defaults(green       = False)
parser.set_defaults(defgrad     = ['f'])

(options,filenames) = parser.parse_args()

stretches = []
stretch = {}
strains = []

if options.right: stretches.append('U')
if options.left:  stretches.append('V')
if options.logarithmic: strains.append('ln')
if options.biot: strains.append('Biot')
if options.green: strains.append('Green')

datainfo = {                                                                                        # list of requested labels per datatype
             'defgrad':     {'len':9,
                             'label':[]},
           }

datainfo['defgrad']['label'] = options.defgrad

# ------------------------------------------ setup file handles ---------------------------------------
files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

# ------------------------------------------ loop over input files ---------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  table = damask.ASCIItable(file['input'],file['output'],False)                                     # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.info_append(string.replace(scriptID,'\n','\\n') + '\t' + ' '.join(sys.argv[1:]))

  active = []
  column = defaultdict(dict)

  for label in datainfo['defgrad']['label']:
    key = '1_%s'%label
    if key not in table.labels:
      sys.stderr.write('column %s not found...\n'%key)
    else:
      active.append(label)
      column[label] = table.labels.index(key)

# ------------------------------------------ assemble header ---------------------------------------

  for label in active:
    for theStretch in stretches:
      for theStrain in strains:
        table.labels_append(['%i_%s(%s)%s'%(i+1,theStrain,theStretch,
                             {True: label,False: ''}[label!='f'])for i in xrange(9)])               # extend ASCII header with new labels  
  table.head_write()

# ------------------------------------------ process data ----------------------------------------  
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    for label in active:                                                                             # loop over all requested norms
      F = np.array(map(float,table.data[column[label]:
                                        column[label]+datainfo['defgrad']['len']]),'d').reshape(3,3)
      (U,S,Vh) = np.linalg.svd(F)
      R = np.dot(U,Vh)
      stretch['U'] = np.dot(np.linalg.inv(R),F)
      stretch['V'] = np.dot(F,np.linalg.inv(R))
      for theStretch in stretches:
        for i in range(9):
          if abs(stretch[theStretch][i%3,i//3]) < 1e-12:                                            # kill nasty noisy data
            stretch[theStretch][i%3,i//3] = 0.0
        (D,V) = np.linalg.eig(stretch[theStretch])                                                  # eigen decomposition (of symmetric matrix)
        for i,eigval in enumerate(D):
          if eigval < 0.0:                                                                          # flip negative eigenvalues
            D[i] = -D[i]
            V[:,i] = -V[:,i]
        if np.dot(V[:,i],V[:,(i+1)%3]) != 0.0:                                                      # check each vector for orthogonality
            V[:,(i+1)%3] = np.cross(V[:,(i+2)%3],V[:,i])                                            # correct next vector
            V[:,(i+1)%3] /= np.sqrt(np.dot(V[:,(i+1)%3],V[:,(i+1)%3].conj()))                       # and renormalize (hyperphobic?)
        for theStrain in strains:
          d = operator(theStretch,theStrain,D)                                                      # operate on eigenvalues of U or V
          eps = (np.dot(V,np.dot(np.diag(d),V.T)).real).reshape(9)                                  # build tensor back from eigenvalue/vector basis

          table.data_append(list(eps))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result ---------------------------------------  
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  file['input'].close()                                                                             # close input ASCII table (works for stdin)
  file['output'].close()                                                                            # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
