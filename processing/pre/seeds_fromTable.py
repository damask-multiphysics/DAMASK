#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,itertools
import numpy as np
from optparse import OptionParser
from collections import defaultdict
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
identifiers = {
        'grid':   ['a','b','c'],
        'size':   ['x','y','z'],
        'origin': ['x','y','z'],
          }
mappings = {
        'grid':            lambda x: int(x),
        'size':            lambda x: float(x),
        'origin':          lambda x: float(x),
        'homogenization':  lambda x: int(x),
        'microstructures': lambda x: int(x),
          }

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Create seed file by taking microstructure indices from given ASCIItable column.
White and black-listing of microstructure indices is possible.

Examples:
--white 1,2,5 --index grainID isolates grainID entries of value 1, 2, and 5;
--black 1 --index grainID takes all grainID entries except for value 1.

""", version = scriptID)

parser.add_option('-p', '--positions',   dest = 'pos', metavar = 'string',
                  help = 'coordinate label')
parser.add_option('--boundingbox',  dest = 'box', type = 'float', nargs = 6, metavar = ' '.join(['float']*6),
                  help = 'min (x,y,z) and max (x,y,z) coordinates of bounding box [auto]')
parser.add_option('-i', '--index',  dest = 'index', type = 'string', metavar = 'string',
                  help = 'microstructure index label')
parser.add_option('-w','--white',   dest = 'whitelist', action = 'extend',
                  help = 'white list of microstructure indices', metavar = '<LIST>')
parser.add_option('-b','--black',   dest = 'blacklist', action = 'extend',
                  help = 'black list of microstructure indices', metavar = '<LIST>')
parser.set_defaults(pos = 'pos',
                    index ='microstructure',
                   )

(options,filenames) = parser.parse_args()

if options.whitelist != None: options.whitelist = map(int,options.whitelist)
if options.blacklist != None: options.blacklist = map(int,options.blacklist)

# --- loop over input files -------------------------------------------------------------------------

if filenames == []:
  filenames = ['STDIN']

for name in filenames:
  if name == 'STDIN':
    file = {'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m\n')
  else:
    if not os.path.exists(name): continue
    file = {'name':name,
            'input':open(name),
            'output':open(os.path.splitext(name)[0]+ \
                         ('' if options.label == None else '_'+options.label)+ \
                         '.png','w'),
            'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  
  table = damask.ASCIItable(file['input'],file['output'],
                            buffered = False)                                                       # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info

# ------------------------------------------ process data ------------------------------------------

  errors = []
  
  missing_labels = table.data_readArray(options.pos,options.label)
  if len(missing_labels) > 0:
    errors.append('column%s %s not found'%('s' if len(missing_labels) > 1 else '',
                                                ', '.join(missing_labels)))
  
  for label, dim in {options.pos: 3,
                     options.label: 1}.iteritems():
    if table.label_dimension(label) != dim:
      errors.append('column %s has wrong dimension'%label)
  
  if errors != []:
    file['croak'].write('\n'.join(errors))
    table.close(dismiss = True)                                                                     # close ASCII table file handles and delete output file
    continue
  
#--- finding bounding box ------------------------------------------------------------------------------------
  boundingBox = np.array((np.amin(table.data[:,0:3],axis = 0),np.amax(table.data[:,0:3],axis = 0)))
  if options.box:
    boundingBox[0,:] = np.minimum(options.box[0:3],boundingBox[0,:])
    boundingBox[1,:] = np.maximum(options.box[3:6],boundingBox[1,:])

#--- rescaling coordinates ------------------------------------------------------------------------------------
  table.data[:,0:3] -= boundingBox[0,:]
  table.data[:,0:3] /= boundingBox[1,:]-boundingBox[0,:]


#--- filtering of grain voxels ------------------------------------------------------------------------------------
  mask = np.logical_and(\
         np.ones_like(table.data[:,3],bool) \
          if options.whitelist == None \
          else              np.in1d(table.data[:,3].ravel(), options.whitelist).reshape(table.data[:,3].shape),
         np.ones_like(table.data[:,3],bool) \
          if options.blacklist == None \
          else    np.invert(np.in1d(table.data[:,3].ravel(), options.blacklist).reshape(table.data[:,3].shape))
          )
  table.data = table.data[mask]

# ------------------------------------------ output result ---------------------------------------  

# ------------------------------------------ assemble header ---------------------------------------  

  table.info = [
                scriptID,
                'size %s'%(' '.join(list(itertools.chain.from_iterable(zip(['x','y','z'],
                                                                       map(str,boundingBox[1,:]-boundingBox[0,:])))))),
               ]
  table.labels_clear()
  table.labels_append(['1_coords','2_coords','3_coords','microstructure'])                          # implicitly switching label processing/writing on
  table.head_write()
  
  table.data_writeArray()
  table.output_flush()
  
  table.close()                                                                                     # close ASCII tables
