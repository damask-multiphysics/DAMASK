#!/usr/bin/env python

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


parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
Create hexahedral voxels around points in an ASCIItable.
""" + string.replace(scriptID,'\n','\\n')
)

parser.add_option('-p', '--positions',   dest='pos', type='string',
                  help = 'coordinate label')
parser.add_option('-s', '--size',   dest='size', type='float', nargs=3,
                  help = 'x,y,z size of voxel')

parser.set_defaults(pos = 'pos')

(options, filenames) = parser.parse_args()

if options.size == None:
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

# --------------- figure out columns to process
  active = {}
  column = {}
  head = []

  for datatype,info in datainfo.items():
    for label in info['label']:
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

  hexPoints = numpy.array([[-1,-1,-1],
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

  minD = numpy.array(options.size,dtype=float)
#   for i in xrange(3):
#     coords = numpy.unique(table.data[:,i])
#     minD[i] = coords[-1]-coords[0]
#     for j in xrange(len(coords)-1):
#       d = coords[j+1]-coords[j]
#       if d < minD[i]:
#         minD[i] = d
    
  for p in table.data:
    for i,h in enumerate(hexPoints):
      id = Points.InsertNextPoint(p+h*minD/2.)
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

  file['input'].close()                                                   # close input ASCII table
