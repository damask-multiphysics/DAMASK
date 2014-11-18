#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,vtk
import numpy as np
import damask
from optparse import OptionParser

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = scriptID.split()[1][:-3]
     
#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
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


parser = OptionParser(option_class=damask.extendableOption, usage='%prog [seedsfile[s]]', description = """
Produce VTK point mesh from seeds file

""", version = scriptID)

parser.add_option('-g', '--grid', dest='grid', type='int', nargs = 3, metavar='int int int', \
                  help='a,b,c grid of hexahedral box [from seeds file]')
parser.add_option('-s', '--size', dest='size', type='float', nargs = 3, metavar='float float float',\
                  help='x,y,z size of hexahedral box [1.0 along largest grid point number]')

parser.set_defaults(grid = [0,0,0])
parser.set_defaults(size = [0.0,0.0,0.0])

(options, filenames) = parser.parse_args()

#--- setup file handles --------------------------------------------------------------------------  
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
                    'output':sys.stdout,
                    'croak':sys.stdout,
                    })

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  theTable = damask.ASCIItable(file['input'],file['output'])
  theTable.head_read()
  
#--- interpret header ----------------------------------------------------------------------------
  info = {
          'grid':   np.zeros(3,'i'),
          'size':   np.zeros(3,'d'),
          'origin': np.zeros(3,'d'),
         }

  for header in theTable.info:
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
                      'origin   x y z:  %s\n'%(' : '.join(map(str,info['origin']))))
  if 0 not in options.grid:                                                                         # user-specified grid
    info['grid'] = np.array(options.grid)
  if np.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
#    continue

  for i in xrange(3):
    if info['size'][i] <= 0.0:                                                                      # any invalid size?
      info['size'][i] = float(info['grid'][i])/max(info['grid'])
      file['croak'].write('rescaling size %s...\n'%{0:'x',1:'y',2:'z'}[i])


#--- read data --------------------------------------------------------------------------------
  theTable.data_readArray(['x','y','z','microstructure'])
  theTable.data[:,0:3] *= info['size']
  theTable.data[:,0:3] += info['origin']

#--- generate grid --------------------------------------------------------------------------------
  grid = vtk.vtkUnstructuredGrid()
  pts = vtk.vtkPoints()

#--- process microstructure information --------------------------------------------------------------
  IDs = vtk.vtkIntArray()
  IDs.SetNumberOfComponents(1)
  IDs.SetName("GrainID")

  for item in theTable.data:
    pid = pts.InsertNextPoint(item[0:3])
    pointIds = vtk.vtkIdList()
    pointIds.InsertId(0, pid)
    grid.InsertNextCell(1, pointIds)
    IDs.InsertNextValue(int(item[3]))

  grid.SetPoints(pts)
  grid.GetCellData().AddArray(IDs)

#--- write data -----------------------------------------------------------------------------------
  if file['name'] == 'STDIN':
    writer = vtk.vtkUnstructuredGridWriter()
    writer.WriteToOutputStringOn()
    writer.SetFileTypeToASCII()
    writer.SetHeader('# powered by '+scriptID)
    if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput(grid)
    else:
      writer.SetInputData(grid)
    writer.Write()
    sys.stdout.write(writer.GetOutputString()[0:writer.GetOutputStringLength()])
  else:
    (head,tail) = os.path.split(file['name'])
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    writer.SetFileName(os.path.join(head,'seeds_'+os.path.splitext(tail)[0]
                                                   +'.'+writer.GetDefaultFileExtension()))
    if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput(grid)
    else:
      writer.SetInputData(grid)
    writer.Write()
