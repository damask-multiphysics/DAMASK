#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,numpy,vtk
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP

#--------------------------------------------------------------------------------------------------
class extendedOption(Option):
#--------------------------------------------------------------------------------------------------
# used for definition of new option parser action 'extend', which enables to take multiple option arguments
# taken from online tutorial http://docs.python.org/library/optparse.html
    
    ACTIONS = Option.ACTIONS + ("extend",)
    STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
    TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
    ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

    def take_action(self, action, dest, opt, value, values, parser):
        if action == "extend":
            lvalue = value.split(",")
            values.ensure_value(dest, []).extend(lvalue)
        else:
            Option.take_action(self, action, dest, opt, value, values, parser)

     
#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
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

parser = OptionParser(option_class=extendedOption, usage='%prog [geomfile[s]]', description = """
Produce VTK rectilinear mesh of structure data from geom description

""" + string.replace('$Id$','\n','\\n')
)

(options, filenames) = parser.parse_args()

#--- setup file handles --------------------------------------------------------------------------  
files = []
if filenames == []:
  files.append({'name':'STDIN',
                'input':sys.stdin,
                'croak':sys.stderr,
               })
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name,
                    'input':open(name),
                    'croak':sys.stdout,
                    })

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write(file['name']+'\n')

  firstline = file['input'].readline()
  m = re.search('(\d+)\s*head', firstline.lower())
  if m:
    headerlines = int(m.group(1))
    headers  = [file['input'].readline() for i in range(headerlines)]
  else:
    headerlines = 1
    headers = firstline

  content = file['input'].readlines()
  file['input'].close()

#--- interprete header ----------------------------------------------------------------------------
  info = {
          'grid':   numpy.zeros(3,'i'),
          'size':   numpy.zeros(3,'d'),
          'origin': numpy.zeros(3,'d'),
          'microstructures': 0,
          'homogenization':  0
         }

  for header in headers:
    headitems = map(str.lower,header.split())
    if headitems[0] == 'resolution': headitems[0] = 'grid'
    if headitems[0] == 'dimension':  headitems[0] = 'size'
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

  if numpy.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    sys.exit()
  if numpy.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    sys.exit()


#--- generate grid --------------------------------------------------------------------------------
  grid = vtk.vtkRectilinearGrid()
  grid.SetDimensions([x+1 for x in info['grid']])
  temp = [] 
  for i in xrange(3):
    temp.append(vtk.vtkDoubleArray())
    temp[i].SetNumberOfTuples(info['grid'][i]+1)
    for j in range(info['grid'][i]+1):
      temp[i].InsertTuple1(j,j*info['size'][i]/info['grid'][i]+info['origin'][i])
      if i == 0: grid.SetXCoordinates(temp[0])
      if i == 1: grid.SetYCoordinates(temp[1])
      if i == 2: grid.SetZCoordinates(temp[2])

#--- read microstructure information --------------------------------------------------------------
  structure = vtk.vtkIntArray()
  structure.SetName('Microstructures')
  for line in content:  
    items = line.split()
    if len(items) > 2:
      if   items[1].lower() == 'of': items = [int(items[2])]*int(items[0])
      elif items[1].lower() == 'to': items = xrange(int(items[0]),1+int(items[2]))
      else:                            items = map(int,items)
    else:                              items = map(int,items)

    for item in items:
      structure.InsertNextValue(item)

  grid.GetCellData().AddArray(structure)

#--- write data -----------------------------------------------------------------------------------
  if file['name'] == 'STDIN':
    outWriter = vtk.vtkRectilinearGridWriter()
    outWriter.WriteToOutputStringOn()
    outWriter.SetFileTypeToASCII()
    outWriter.SetHeader('# powered by $Id$')
    outWriter.SetInput(grid)
    outWriter.Write()
    sys.stdout.write(outWriter.GetOutputString()[0:outWriter.GetOutputStringLength()])
  else:
    (head,tail) = os.path.split(file['name'])
    outWriter = vtk.vtkXMLRectilinearGridWriter()
    outWriter.SetDataModeToBinary()
    outWriter.SetCompressorTypeToZLib()
    outWriter.SetFileName(os.path.join(head,'mesh_'+os.path.splitext(tail)[0]
                                                   +'.'+outWriter.GetDefaultFileExtension()))
    outWriter.SetInput(grid)
    outWriter.Write()
