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

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Create seed file by taking microstructure indices from given ASCIItable column.
White and black-listing of microstructure indices is possible.

Examples:
--white 1,2,5 --index grainID isolates grainID entries of value 1, 2, and 5;
--black 1 --index grainID takes all grainID entries except for value 1.

""", version = scriptID)

parser.add_option('-p', '--positions',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'coordinate label [%default]')
parser.add_option('--boundingbox',
                  dest = 'box',
                  type = 'float', nargs = 6, metavar = ' '.join(['float']*6),
                  help = 'min (x,y,z) and max (x,y,z) coordinates of bounding box [tight]')
parser.add_option('-i', '--index',
                  dest = 'index',
                  type = 'string', metavar = 'string',
                  help = 'microstructure index label [%default]')
parser.add_option('-w','--white',
                  dest = 'whitelist',
                  action = 'extend', metavar = '<int LIST>',
                  help = 'whitelist of microstructure indices')
parser.add_option('-b','--black',
                  dest = 'blacklist',
                  action = 'extend', metavar = '<int LIST>',
                  help = 'blacklist of microstructure indices')

parser.set_defaults(pos = 'pos',
                    index ='microstructure',
                   )

(options,filenames) = parser.parse_args()

if options.whitelist != None: options.whitelist = map(int,options.whitelist)
if options.blacklist != None: options.blacklist = map(int,options.blacklist)

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = ['STDIN']

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              outname = os.path.splitext(name)[0]+'.seeds' if name else name,
                              buffered = False)
  except: continue
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name else ''))

  table.head_read()                                                                                 # read ASCII header info

# ------------------------------------------ sanity checks ---------------------------------------  

  missing_labels = table.data_readArray([options.pos,options.index])

  errors = []
  if len(missing_labels) > 0:
    errors.append('column{} {} not found'.format('s' if len(missing_labels) > 1 else '',
                                                ', '.join(missing_labels)))
  for label, dim in {options.pos: 3,
                     options.index: 1}.iteritems():
    if table.label_dimension(label) != dim:
      errors.append('column {} has wrong dimension'.format(label))
  
  if errors != []:
    table.croak(errors)
    table.close(dismiss = True)                                                                     # close ASCII table file handles and delete output file
    continue

# ------------------------------------------ process data ------------------------------------------
  
# --- finding bounding box -------------------------------------------------------------------------

  boundingBox = np.array((np.amin(table.data[:,0:3],axis = 0),np.amax(table.data[:,0:3],axis = 0)))
  if options.box:
    boundingBox[0,:] = np.minimum(options.box[0:3],boundingBox[0,:])
    boundingBox[1,:] = np.maximum(options.box[3:6],boundingBox[1,:])

# --- rescaling coordinates ------------------------------------------------------------------------

  table.data[:,0:3] -= boundingBox[0,:]
  table.data[:,0:3] /= boundingBox[1,:]-boundingBox[0,:]

# --- filtering of grain voxels --------------------------------------------------------------------

  mask = np.logical_and(\
         np.ones_like(table.data[:,3],bool) \
          if options.whitelist == None \
          else              np.in1d(table.data[:,3].ravel(), options.whitelist).reshape(table.data[:,3].shape),
         np.ones_like(table.data[:,3],bool) \
          if options.blacklist == None \
          else    np.invert(np.in1d(table.data[:,3].ravel(), options.blacklist).reshape(table.data[:,3].shape))
          )
  table.data = table.data[mask]

# ------------------------------------------ assemble header ---------------------------------------  

  table.info = [
                scriptID,
                'size %s'%(' '.join(list(itertools.chain.from_iterable(zip(['x','y','z'],
                                                                       map(str,boundingBox[1,:]-boundingBox[0,:])))))),
               ]
  table.labels_clear()
  table.labels_append(['1_pos','2_pos','3_pos','microstructure'])                                   # implicitly switching label processing/writing on
  table.head_write()
  
# ------------------------------------------ output result ---------------------------------------  

  table.data_writeArray()
  table.close()                                                                                     # close ASCII tables
