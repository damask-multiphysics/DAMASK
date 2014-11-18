#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,vtk
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]
  
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
        'homogenization':  lambda x: int(x),
        'microstructures': lambda x: int(x),
          }

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [file[s]]', description = """
Produce VTK rectilinear mesh of structure data from geom description

""", version = scriptID)

parser.add_option('-n','--nodata',      dest='data', action='store_false',
                                        help='omit microstructure data, just generate mesh')

parser.set_defaults(data = True)

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
  file['croak'].write('\033[1m' + scriptName + '\033[0m: ' + (file['name'] if file['name'] != 'STDIN' else '') + '\n')

  theTable = damask.ASCIItable(file['input'],file['output'],labels=False)
  theTable.head_read()

#--- interpret header ----------------------------------------------------------------------------
  info = {
          'grid':   np.zeros(3,'i'),
          'size':   np.zeros(3,'d'),
          'origin': np.zeros(3,'d'),
          'homogenization':  0,
          'microstructures': 0,
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
                      'origin   x y z:  %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization:  %i\n'%info['homogenization'] + \
                      'microstructures: %i\n'%info['microstructures'])

  if np.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    continue
  if np.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    continue


#--- generate grid --------------------------------------------------------------------------------
  grid = vtk.vtkRectilinearGrid()
  grid.SetDimensions([x+1 for x in info['grid']])
  for i in xrange(3):
    temp = vtk.vtkDoubleArray()
    temp.SetNumberOfTuples(info['grid'][i]+1)
    for j in xrange(info['grid'][i]+1):
      temp.InsertTuple1(j,j*info['size'][i]/info['grid'][i]+info['origin'][i])
    if   i == 0: grid.SetXCoordinates(temp)
    elif i == 1: grid.SetYCoordinates(temp)
    elif i == 2: grid.SetZCoordinates(temp)

#--- read microstructure information --------------------------------------------------------------

  if options.data:
    structure = vtk.vtkIntArray()
    structure.SetName('Microstructures')

    while theTable.data_read():
      items = theTable.data
      if len(items) > 2:
        if   items[1].lower() == 'of': items = [int(items[2])]*int(items[0])
        elif items[1].lower() == 'to': items = xrange(int(items[0]),1+int(items[2]))
        else:                          items = map(int,items)
      else:                            items = map(int,items)

      for item in items:
        structure.InsertNextValue(item)

    grid.GetCellData().AddArray(structure)

#--- write data -----------------------------------------------------------------------------------
  if file['name'] == 'STDIN':
    writer = vtk.vtkRectilinearGridWriter()
    writer.WriteToOutputStringOn()
    writer.SetFileTypeToASCII()
    writer.SetHeader('# powered by '+scriptID)
    if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput(grid)
    else:
      writer.SetInputData(grid)
    writer.Write()
    file['output'].write(writer.GetOutputString()[0:writer.GetOutputStringLength()])
  else:
    (dir,file) = os.path.split(file['name'])
    writer = vtk.vtkXMLRectilinearGridWriter()
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    writer.SetFileName(os.path.join(dir,'mesh_'+os.path.splitext(file)[0]
                                                   +'.'+writer.GetDefaultFileExtension()))
    if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput(grid)
    else:
      writer.SetInputData(grid)
    writer.Write()
