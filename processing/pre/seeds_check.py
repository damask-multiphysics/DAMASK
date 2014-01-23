#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,numpy,vtk
import damask
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP

scriptID = '$Id$'
scriptName = scriptID.split()[1]

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


parser = OptionParser(option_class=extendedOption, usage='%prog [seedsfile[s]]', description = """
Produce VTK point mesh from seeds file

""" + string.replace(scriptID,'\n','\\n')
)

parser.add_option('-g', '--grid', dest='grid', type='int', nargs = 3, \
                  help='a,b,c grid of hexahedral box [from seeds file]')
parser.add_option('-s', '--size', dest='size', type='float', nargs = 3, \
                  help='x,y,z size of hexahedral box [1.0 along largest grid point number]')

parser.set_defaults(grid = [0,0,0])
parser.set_defaults(size  = [0.0,0.0,0.0])

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
  
  coords = theTable.data_asArray(['x','y','z'])
  if theTable.labels_index('microstructure') != -1:
    grain = theTable.data_asArray(['microstructure'])
    grainIDs = numpy.unique(grain).astype('i')
  else:
    grain = 1+numpy.arange(len(coords))
    grainIDs = grain

#--- interpret header ----------------------------------------------------------------------------
  info = {
          'grid':   numpy.zeros(3,'i'),
          'size':   numpy.zeros(3,'d'),
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
                      'size     x y z:  %s\n'%(' x '.join(map(str,info['size']))))

  if 0 not in options.grid:                                                                         # user-specified grid
    info['grid'] = numpy.array(options.grid)
  if numpy.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    continue

  for i in xrange(3):
    if info['size'][i] <= 0.0:                                                                      # any invalid size?
      info['size'][i] = float(info['grid'][i])/max(info['grid'])
      file['croak'].write('rescaling size %s...\n'%{0:'x',1:'y',2:'z'}[i])


#--- generate grid --------------------------------------------------------------------------------
  grid = vtk.vtkUnstructuredGrid()
  pts = vtk.vtkPoints()

#--- process microstructure information --------------------------------------------------------------
  IDs = vtk.vtkIntArray()
  IDs.SetNumberOfComponents(1)
  IDs.SetName("GrainID")

  for coord,ID in zip(coords,grainIDs):
    pid = pts.InsertNextPoint(coord)
    pointIds = vtk.vtkIdList()
    pointIds.InsertId(0, pid)
    grid.InsertNextCell(1, pointIds)
    IDs.InsertNextValue(ID)

  grid.SetPoints(pts)
  grid.GetCellData().AddArray(IDs)

#--- write data -----------------------------------------------------------------------------------
  if file['name'] == 'STDIN':
    outWriter = vtk.vtkUnstructuredGridWriter()
    outWriter.WriteToOutputStringOn()
    outWriter.SetFileTypeToASCII()
    outWriter.SetHeader('# powered by '+scriptID)
    outWriter.SetInputData(grid)
    outWriter.Write()
    sys.stdout.write(outWriter.GetOutputString()[0:outWriter.GetOutputStringLength()])
  else:
    (head,tail) = os.path.split(file['name'])
    outWriter = vtk.vtkXMLUnstructuredGridWriter()
    outWriter.SetDataModeToBinary()
    outWriter.SetCompressorTypeToZLib()
    outWriter.SetFileName(os.path.join(head,'seeds_'+os.path.splitext(tail)[0]
                                                   +'.'+outWriter.GetDefaultFileExtension()))
    outWriter.SetInputData(grid)
    outWriter.Write()
