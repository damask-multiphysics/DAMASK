#!/usr/bin/env python3

import os
import argparse

import numpy as np

import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
parser = argparse.ArgumentParser()

#ToDo:  We need to decide on a way of handling arguments of variable lentght
#https://stackoverflow.com/questions/15459997/passing-integer-lists-to-python

#parser.add_argument('--version', action='version', version='%(prog)s {}'.format(scriptID))
parser.add_argument('filenames', nargs='+',
                    help='DADF5 files')
parser.add_argument('-d','--dir', dest='dir',default='postProc',metavar='string',
                    help='name of subdirectory to hold output')
parser.add_argument('--mat', nargs='+',
                    help='labels for materialpoint/homogenization',dest='mat')
parser.add_argument('--con', nargs='+',
                    help='labels for constituent/crystallite/constitutive',dest='con')

options = parser.parse_args()

if options.mat is None: options.mat=[]
if options.con is None: options.con=[]

# --- loop over input files ------------------------------------------------------------------------

for filename in options.filenames:
  results = damask.DADF5(filename)
  
  if not results.structured: continue
  delta = results.size/results.grid*0.5
  x, y, z = np.meshgrid(np.linspace(delta[2],results.size[2]-delta[2],results.grid[2]),
                        np.linspace(delta[1],results.size[1]-delta[1],results.grid[1]),
                        np.linspace(delta[0],results.size[0]-delta[0],results.grid[0]),
                        indexing = 'ij')

  coords = np.concatenate((z[:,:,:,None],y[:,:,:,None],x[:,:,:,None]),axis = 3) 
  
  for i,inc in enumerate(results.increments):
    print('Output step {}/{}'.format(i+1,len(results.increments)))

    header = '1 header\n'
    
    data = np.array([inc['inc'] for j in range(np.product(results.grid))]).reshape([np.product(results.grid),1])
    header+= 'inc'

    coords = coords.reshape([np.product(results.grid),3])
    data = np.concatenate((data,coords),1)
    header+=' 1_pos 2_pos 3_pos'
        
    results.visible['increments'] = [inc]

    for label in options.con:
      for o in results.constituent_output_iter():
        for c in results.constituent_iter():
          x = results.get_dataset_location(label)
          if len(x) == 0:
            continue
          label = x[0].split('/')[-1]
          array = results.read_dataset(x,0)
          d = int(np.product(np.shape(array)[1:]))
          array = np.reshape(array,[np.product(results.grid),d])
          data = np.concatenate((data,array),1)
          
          if d>1:
            header+= ''.join([' {}_{}'.format(j+1,label) for j in range(d)])
          else:
            header+=' '+label
            
    for label in options.mat:
      for o in results.materialpoint_output_iter():
        for m in results.materialpoint_iter():
          x = results.get_dataset_location(label)
          if len(x) == 0:
            continue
          label = x[0].split('/')[-1]
          array = results.read_dataset(x,0)
          d = int(np.product(np.shape(array)[1:]))
          array = np.reshape(array,[np.product(results.grid),d])
          data = np.concatenate((data,array),1)
          
          if d>1:
            header+= ''.join([' {}_{}'.format(j+1,label) for j in range(d)])
          else:
            header+=' '+label

    dirname  = os.path.abspath(os.path.join(os.path.dirname(filename),options.dir))
    try:
      os.mkdir(dirname)
    except FileExistsError:
      pass
    file_out = '{}_inc{:04d}.txt'.format(filename.split('.')[0],inc['inc'])
    np.savetxt(os.path.join(dirname,file_out),data,header=header,comments='')
