#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math,string
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

oversampling = 2.

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
Positions a geometric object within the (three-dimensional) canvas of a spectral geometry description.
Depending on the sign of the dimension parameters, these objects can be boxes, cylinders, or ellipsoids.

""", version = scriptID)

parser.add_option('-c', '--center',     dest='center', type='int', nargs = 3, metavar=' '.join(['int']*3),
                  help='a,b,c origin of primitive %default')
parser.add_option('-d', '--dimension',  dest='dimension', type='int', nargs = 3, metavar=' '.join(['int']*3),
                  help='a,b,c extension of hexahedral box; negative values are diameters')
parser.add_option('-f', '--fill',       dest='fill', type='int', metavar = 'int',
                  help='grain index to fill primitive. "0" selects maximum microstructure index + 1 [%default]')
parser.add_option('-q', '--quaternion', dest='quaternion', type='float', nargs = 4, metavar=' '.join(['float']*4),
                  help = 'rotation of primitive as quaternion')
parser.add_option('-a', '--angleaxis',  dest='angleaxis', nargs = 4, metavar=' '.join(['float']*4),
                  help = 'rotation of primitive as angle and axis')
parser.add_option(     '--degrees',     dest='degrees', action='store_true',
                  help = 'angle is given in degrees [%default]')

parser.set_defaults(center = [0,0,0],
                    fill = 0,
                    quaternion = [],
                    angleaxis = [],
                    degrees = False,
                   )

(options, filenames) = parser.parse_args()

if options.angleaxis != []:
  options.angleaxis = map(float,options.angleaxis)
  rotation = damask.Quaternion().fromAngleAxis(np.radians(options.angleaxis[0]) if options.degrees else options.angleaxis[0],
                                               options.angleaxis[1:4]).conjugated()
elif options.quaternion != []:
  options.quaternion = map(float,options.quaternion)
  rotation = damask.Quaternion(options.quaternion).conjugated()
else:
  rotation = damask.Quaternion().conjugated()

options.center = np.array(options.center)
invRotation = rotation.conjugated()                                                             # rotation of gridpos into primitive coordinate system

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

#--- read data ------------------------------------------------------------------------------------

  microstructure = table.microstructure_read(info['grid'])                                          # read microstructure

# --- do work ------------------------------------------------------------------------------------

  newInfo = {
             'microstructures': 0,
            }


  if options.fill == 0:
    options.fill = microstructure.max()+1

  microstructure = microstructure.reshape(info['grid'],order='F')

  if options.dimension != None:
    mask = (np.array(options.dimension) < 0).astype(float)                                       # zero where positive dimension, otherwise one
    dim = abs(np.array(options.dimension))                                                       # dimensions of primitive body
    pos = np.zeros(3,dtype='float')
#    hiresPrimitive = np.zeros((2*dim[0],2*dim[1],2*dim[2],3))                                   # primitive discretized at twice the grid resolution
    for     i,pos[0] in enumerate(np.arange(-dim[0]/oversampling,(dim[0]+1)/oversampling,1./oversampling)):
      for   j,pos[1] in enumerate(np.arange(-dim[1]/oversampling,(dim[1]+1)/oversampling,1./oversampling)):
        for k,pos[2] in enumerate(np.arange(-dim[2]/oversampling,(dim[2]+1)/oversampling,1./oversampling)):
          gridpos = np.floor(rotation*pos)                                                       # rotate and lock into spacial grid
          primPos = invRotation*gridpos                                                             # rotate back to primitive coordinate system
          if np.dot(mask*primPos/dim,mask*primPos/dim) <= 0.25 and \
             np.all(abs((1.-mask)*primPos/dim) <= 0.5):                                          # inside ellipsoid and inside box
             microstructure[(gridpos[0]+options.center[0])%info['grid'][0],
                            (gridpos[1]+options.center[1])%info['grid'][1],
                            (gridpos[2]+options.center[2])%info['grid'][2]] = options.fill          # assign microstructure index

  newInfo['microstructures'] = microstructure.max()


# --- report ---------------------------------------------------------------------------------------

  remarks = []
  if (    newInfo['microstructures'] != info['microstructures']): remarks.append('--> microstructures: %i'%newInfo['microstructures'])
  if remarks != []: damask.util.croak(remarks)

#--- write header ---------------------------------------------------------------------------------

  table.info_clear()
  table.info_append([
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=info['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=info['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=info['origin']),
    "homogenization\t{homog}".format(homog=info['homogenization']),
    "microstructures\t{microstructures}".format(microstructures=newInfo['microstructures']),
    extra_header
    ])
  table.labels_clear()
  table.head_write()
  table.output_flush()
    
# --- write microstructure information ------------------------------------------------------------

  formatwidth = int(math.floor(math.log10(microstructure.max())+1))
  table.data = microstructure.reshape((info['grid'][0],info['grid'][1]*info['grid'][2]),order='F').transpose()
  table.data_writeArray('%%%ii'%(formatwidth),delimiter = ' ')

#--- output finalization --------------------------------------------------------------------------

  table.close()
