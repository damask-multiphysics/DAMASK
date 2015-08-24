#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,vtk
import damask
from optparse import OptionParser

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Produce a VTK point cloud dataset based on coordinates given in an ASCIItable.

""", version = scriptID)

parser.add_option('-p','--position',
                  dest = 'position',
                  type = 'string', metavar = 'string',
                  help = 'column label for coordinates [%default]')

parser.set_defaults(position = 'pos',
                   )

(options, filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False, readonly = True)
  except: continue
  table.croak(damask.util.emph(scriptName)+(': '+name if name else ''))

# --- interpret header ----------------------------------------------------------------------------

  table.head_read()

  errors =  []
  remarks = []
  if table.label_dimension(options.position) != 3: errors.append('columns "{}" have dimension {}'.format(options.position,
                                                                                                         table.label_dimension(options.position)))
  if remarks != []: table.croak(remarks)
  if errors  != []:
    table.croak(errors)
    table.close(dismiss=True)
    continue

# ------------------------------------------ process data ---------------------------------------  

  table.data_readArray(options.position)

  Polydata = vtk.vtkPolyData()
  Points = vtk.vtkPoints()
  Vertices = vtk.vtkCellArray()

  for p in table.data:
    id = Points.InsertNextPoint(p)
    Vertices.InsertNextCell(1)
    Vertices.InsertCellPoint(id)

  Polydata.SetPoints(Points)
  Polydata.SetVerts(Vertices)
  Polydata.Modified()
  if vtk.VTK_MAJOR_VERSION <= 5: Polydata.Update()
 
# ------------------------------------------ output result ---------------------------------------  

  if name:
    (dir,filename) = os.path.split(name)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    writer.SetFileName(os.path.join(dir,os.path.splitext(filename)[0]+'_{}'.format(options.position) \
                                        +'.'+writer.GetDefaultFileExtension()))
    if vtk.VTK_MAJOR_VERSION <= 5: writer.SetInput(Polydata)
    else:                          writer.SetInputData(Polydata)
    writer.Write()

  else:
    writer = vtk.vtkPolyDataWriter()
    writer.WriteToOutputStringOn()
    writer.SetFileTypeToASCII()
    writer.SetHeader('# powered by '+scriptID)
    if vtk.VTK_MAJOR_VERSION <= 5: writer.SetInput(Polydata)
    else:                          writer.SetInputData(Polydata)
    writer.Write()
    sys.stdout.write(writer.GetOutputString()[0:writer.GetOutputStringLength()])

  table.close()                                                     # close input ASCII table
