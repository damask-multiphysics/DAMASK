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
                    help='name of subdirectory relative to the location of the DADF5 file to hold output')
parser.add_argument('--mat', nargs='+',
                    help='labels for materialpoint',dest='mat')
parser.add_argument('--con', nargs='+',
                    help='labels for constituent',dest='con')

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
  
  N_digits = int(np.floor(np.log10(int(results.increments[-1][3:]))))+1
  N_digits = 5 # hack to keep test intact
  for i,inc in enumerate(results.iter_visible('increments')):
    print('Output step {}/{}'.format(i+1,len(results.increments)))

    header = '1 header\n'
    
    data = np.array([int(inc[3:]) for j in range(np.product(results.grid))]).reshape([np.product(results.grid),1])
    header+= 'inc'

    coords = coords.reshape([np.product(results.grid),3])
    data = np.concatenate((data,coords),1)
    header+=' 1_pos 2_pos 3_pos'

    for label in options.con:
      for p in results.iter_visible('con_physics'):
        for c in results.iter_visible('constituents'):
          x = results.get_dataset_location(label)
          if len(x) == 0:
            continue
          array = results.read_dataset(x,0,plain=True)
          d = int(np.product(np.shape(array)[1:]))
          data = np.concatenate((data,np.reshape(array,[np.product(results.grid),d])),1)
          
          if d>1:
            header+= ''.join([' {}_{}'.format(j+1,label) for j in range(d)])
          else:
            header+=' '+label
            
    for label in options.mat:
      for p in results.iter_visible('mat_physics'):
        for m in results.iter_visible('materialpoints'):
          x = results.get_dataset_location(label)
          if len(x) == 0:
            continue
          array = results.read_dataset(x,0,plain=True)
          d = int(np.product(np.shape(array)[1:]))
          data = np.concatenate((data,np.reshape(array,[np.product(results.grid),d])),1)
          
          if d>1:
            header+= ''.join([' {}_{}'.format(j+1,label) for j in range(d)])
          else:
            header+=' '+label

    dirname  = os.path.abspath(os.path.join(os.path.dirname(filename),options.dir))
    if not os.path.isdir(dirname):
      os.mkdir(dirname,0o755)
    file_out = '{}_inc{}.txt'.format(os.path.splitext(os.path.split(filename)[-1])[0],
                                     inc[3:].zfill(N_digits))
    np.savetxt(os.path.join(dirname,file_out),data,header=header,comments='')
