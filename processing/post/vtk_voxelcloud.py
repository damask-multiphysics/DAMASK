#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,vtk
import numpy as np
from optparse import OptionParser
import damask

scriptID = '$Id$'
scriptName = os.path.splitext(scriptID.split()[1])[0]

synonyms = {
        'grid':   ['resolution'],
        'size':   ['dimension'],
          }
identifiers = {
        'grid':    ['a','b','c'],
        'size':    ['x','y','z'],
        'origin':  ['x','y','z'],
          }
mappings = {
        'grid':            lambda x: int(x),
        'size':            lambda x: float(x),
        'origin':          lambda x: float(x),
        'microstructures': lambda x: int(x),
          }

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Create hexahedral voxels around points in an ASCIItable.

""", version = scriptID)

parser.add_option('-p', '--positions',   dest='pos', type='string',
                  help = 'coordinate label')
parser.add_option('-s', '--size',   dest='size', type='float', nargs=3,
                  help = 'x,y,z size of voxel')
parser.add_option('-o', '--origin',   dest='origin', type='float', nargs=3,
                  help = 'x,y,z origin of coordinate system')
parser.add_option('-g', '--geom',   dest='geom', action='store_true',
                  help = 'derive geometry from geom-file header information')

parser.set_defaults(pos = 'pos')
parser.set_defaults(origin = [0.0,0.0,0.0])
parser.set_defaults(geom = False)

(options, filenames) = parser.parse_args()

if options.size == None and not options.geom:
  parser.error('no size sprecified.')

datainfo = {                                                               # list of requested labels per datatype
             'vector':     {'len':3,
                            'label':[]},
           }

if options.pos != None:  datainfo['vector']['label'] += [options.pos]

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
for name in filenames:
  if os.path.exists(name):
    files.append({'name':name, 'input':open(name), 'output':os.path.splitext(name)[0]+'.vtu', 'croak':sys.stderr})

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  table = damask.ASCIItable(file['input'],file['croak'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info
#--- interpret header ----------------------------------------------------------------------------
  info = {
          'grid':   np.zeros(3,'i'),
          'size':   np.zeros(3,'d'),
          'origin': np.zeros(3,'d'),
          'homogenization':  0,
          'microstructures': 0,
         }

  if options.geom:
    for header in table.info:
      headitems = map(str.lower,header.split())
      if len(headitems) == 0: continue
      for synonym,alternatives in synonyms.iteritems():
        if headitems[0] in alternatives: headitems[0] = synonym
      if headitems[0] in mappings.keys():
        if headitems[0] in identifiers.keys():
          for i in xrange(len(identifiers[headitems[0]])):
            info[headitems[0]][i] = \
              mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])
        else:
          info[headitems[0]] = mappings[headitems[0]](headitems[1])

    file['croak'].write('grid     a b c:  %s\n'%(' x '.join(map(str,info['grid']))) + \
                        'size     x y z:  %s\n'%(' x '.join(map(str,info['size']))) + \
                        'origin   x y z:  %s\n'%(' : '.join(map(str,info['origin']))) + \
                        'homogenization:  %i\n'%info['homogenization'] + \
                        'microstructures: %i\n'%info['microstructures'])

    if np.any(info['grid'] < 1):
      file['croak'].write('invalid grid a b c.\n')
      continue
    if np.any(info['size'] <= 0.0):
      file['croak'].write('invalid size x y z.\n')
      continue

  else:
    info['size'] = np.ones(3)
    info['grid'] = info['size'] / options.size
    info['origin'] = options.origin

# --------------- figure out columns to process
  active = {}
  column = {}
  head = []

  for datatype,infos in datainfo.items():
    for label in infos['label']:
      foundIt = False
      for key in ['1_'+label,label]:
        if key in table.labels:
          foundIt = True
          if datatype not in active: active[datatype] = []
          if datatype not in column: column[datatype] = {}
          active[datatype].append(label)
          column[datatype][label] = table.labels.index(key)                 # remember columns of requested data
      if not foundIt:
        file['croak'].write('column %s not found...\n'%label)
        break


# ------------------------------------------ process data ---------------------------------------  

  hexPoints = np.array([[-1,-1,-1],
                           [ 1,-1,-1],
                           [ 1, 1,-1],
                           [-1, 1,-1],
                           [-1,-1, 1],
                           [ 1,-1, 1],
                           [ 1, 1, 1],
                           [-1, 1, 1],
                          ])
  Points = vtk.vtkPoints()
  Hex = vtk.vtkHexahedron()
  uGrid = vtk.vtkUnstructuredGrid()

  table.data_readArray(range(column['vector'][options.pos],\
                             column['vector'][options.pos]+datainfo['vector']['len']))

  table.data[:,0:3] *= info['size']
  table.data[:,0:3] += info['origin']

#  minD = np.array(options.size,dtype=float)
#   for i in xrange(3):
#     coords = np.unique(table.data[:,i])
#     minD[i] = coords[-1]-coords[0]
#     for j in xrange(len(coords)-1):
#       d = coords[j+1]-coords[j]
#       if d < minD[i]:
#         minD[i] = d
    
  for p in table.data:
    for i,h in enumerate(hexPoints):
      id = Points.InsertNextPoint(p+h*info['size']/info['grid']/2.)
      Hex.GetPointIds().SetId(i,id)

    uGrid.InsertNextCell(Hex.GetCellType(), Hex.GetPointIds())

  uGrid.SetPoints(Points)

# ------------------------------------------ output result ---------------------------------------  

  writer = vtk.vtkXMLUnstructuredGridWriter()
  writer.SetDataModeToBinary()
  writer.SetCompressorTypeToZLib()
  writer.SetFileName(file['output'])
  if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput(uGrid)
  else:
      writer.SetInputData(uGrid)
  writer.Write()

  table.input_close()                                                     # close input ASCII table
