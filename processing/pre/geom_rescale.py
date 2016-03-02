#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Scales a geometry description independently in x, y, and z direction in terms of grid and/or size.
Either absolute values or relative factors (like "0.25x") can be used.

""", version = scriptID)

parser.add_option('-g', '--grid',
                  dest = 'grid',
                  type = 'string', nargs = 3, metavar = 'string string string',
                  help = 'a,b,c grid of hexahedral box [unchanged]')
parser.add_option('-s', '--size',
                  dest = 'size',
                  type = 'string', nargs = 3, metavar = 'string string string',
                  help = 'x,y,z size of hexahedral box [unchanged]')
parser.add_option('-r', '--renumber',
                  dest = 'renumber',
                  action = 'store_true',
                  help = 'renumber microstructure indices from 1..N [%default]')

parser.set_defaults(renumber = False,
                    grid = ['0','0','0'],
                    size  = ['0.0','0.0','0.0'],
                   )

(options, filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False, labeled = False)
  except: continue
  damask.util.report(scriptName,name)

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()
  info,extra_header = table.head_getGeom()
  
  damask.util.croak(['grid     a b c:  %s'%(' x '.join(map(str,info['grid']))),
               'size     x y z:  %s'%(' x '.join(map(str,info['size']))),
               'origin   x y z:  %s'%(' : '.join(map(str,info['origin']))),
               'homogenization:  %i'%info['homogenization'],
               'microstructures: %i'%info['microstructures'],
              ])

  errors = []
  if np.any(info['grid'] < 1):    errors.append('invalid grid a b c.')
  if np.any(info['size'] <= 0.0): errors.append('invalid size x y z.')
  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --- read data ------------------------------------------------------------------------------------

  microstructure = table.microstructure_read(info['grid'])                                          # read microstructure

# --- do work ------------------------------------------------------------------------------------
 
  newInfo = {
             'grid':    np.zeros(3,'i'),
             'origin':  np.zeros(3,'d'),
             'microstructures': 0,
            }

  newInfo['grid'] = np.array([{True:round(o*float(n.translate(None,'xX'))), 
                               False: round(float(n.translate(None,'xX')))}[n[-1].lower() == 'x']
                                                  for o,n in zip(info['grid'],options.grid)],'i')
  newInfo['size'] = np.array([{True:      o*float(n.translate(None,'xX')) ,
                               False:       float(n.translate(None,'xX')) }[n[-1].lower() == 'x']
                                                  for o,n in zip(info['size'],options.size)],'d')
  newInfo['grid'] = np.where(newInfo['grid'] <= 0  , info['grid'],newInfo['grid'])
  newInfo['size'] = np.where(newInfo['size'] <= 0.0, info['size'],newInfo['size'])

  multiplicity = []
  for j in xrange(3):
    multiplicity.append([])
    last = 0
    for i in xrange(info['grid'][j]):
      this = int((i+1)*float(newInfo['grid'][j])/info['grid'][j])
      multiplicity[j].append(this-last)
      last = this

  microstructure = microstructure.reshape(info['grid'],order='F')
  microstructure = np.repeat(np.repeat(np.repeat(microstructure,
                   multiplicity[0], axis=0),multiplicity[1], axis=1),multiplicity[2], axis=2)
# --- renumber to sequence 1...Ngrains if requested ------------------------------------------------
#  http://stackoverflow.com/questions/10741346/np-frequency-counts-for-unique-values-in-an-array  

  if options.renumber:
    newID = 0
    for microstructureID,count in enumerate(np.bincount(microstructure.reshape(newInfo['grid'].prod()))):
      if count != 0:
        newID += 1
        microstructure = np.where(microstructure == microstructureID, newID,microstructure).reshape(microstructure.shape)

  newInfo['microstructures'] = microstructure.max()

# --- report ---------------------------------------------------------------------------------------

  remarks = []
  errors  = []

  if (any(newInfo['grid']            != info['grid'])):
    remarks.append('--> grid     a b c:  %s'%(' x '.join(map(str,newInfo['grid']))))
  if (any(newInfo['size']            != info['size'])):
    remarks.append('--> size     x y z:  %s'%(' x '.join(map(str,newInfo['size']))))
  if (    newInfo['microstructures'] != info['microstructures']):
    remarks.append('--> microstructures: %i'%newInfo['microstructures'])

  if np.any(newInfo['grid'] < 1):    errors.append('invalid new grid a b c.')
  if np.any(newInfo['size'] <= 0.0): errors.append('invalid new size x y z.')

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --- write header ---------------------------------------------------------------------------------

  table.info_clear()
  table.info_append([
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=newInfo['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=newInfo['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=info['origin']),
    "homogenization\t{homog}".format(homog=info['homogenization']),
    "microstructures\t{microstructures}".format(microstructures=newInfo['microstructures']),
    extra_header
    ])
  table.labels_clear()
  table.head_write()
  
# --- write microstructure information ------------------------------------------------------------

  formatwidth = int(math.floor(math.log10(microstructure.max())+1))
  table.data = microstructure.reshape((newInfo['grid'][0],newInfo['grid'][1]*newInfo['grid'][2]),order='F').transpose()
  table.data_writeArray('%%%ii'%(formatwidth),delimiter = ' ')
    
# --- output finalization --------------------------------------------------------------------------

  table.close()                                                                                     # close ASCII table
