#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,sys,string,re,numpy,scipy.ndimage,scipy.signal,vtk
import damask
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP

scriptID = '$Id$'
scriptName = scriptID.split()[1]

scalingFactor = { \
                 'm':  {
                        'm':  1e0,
                        'mm': 1e-3,
                        'µm': 1e-6,
                       },
                 'mm':  {
                        'm':  1e+3,
                        'mm': 1e0,
                        'µm': 1e-3,
                       },
                 'µm':  {
                        'm':  1e+6,
                        'mm': 1e+3,
                        'µm': 1e0,
                       },
                }

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
Produce VTK rectilinear grid from Gwyddion dataset exported as text.
""" + string.replace(scriptID,'\n','\\n')
)

parser.add_option('-s', '--scaling',   dest='scaling', type='float',
                  help = 'scaling factor for elevation data [auto]')

parser.set_defaults(scaling = 0.0)

(options, filenames) = parser.parse_args()


# ------------------------------------------ read Gwyddion data ---------------------------------------  

for file in filenames:
  with open(file,'r') as f:
    for line in f:
      pieces = line.split()
      if pieces[0] != '#': break
      if len(pieces) < 2: continue
      if pieces[1] == 'Width:':
        width  = float(pieces[2])
        lateralunit = pieces[3]
      if pieces[1] == 'Height:':
        height = float(pieces[2])
        lateralunit = pieces[3]
      if pieces[1] == 'Value' and pieces[2] == 'units:':
        elevationunit = pieces[3]

    if options.scaling == 0.0:
      options.scaling = scalingFactor[lateralunit][elevationunit]
  
    elevation = numpy.loadtxt(file)*options.scaling
    
    grid = vtk.vtkRectilinearGrid()
    grid.SetDimensions(elevation.shape[1],elevation.shape[0],1)

    xCoords = vtk.vtkDoubleArray()
    for x in numpy.arange(0.0,width,width/elevation.shape[1],'d'):
      xCoords.InsertNextValue(x)
    yCoords = vtk.vtkDoubleArray()
    for y in numpy.arange(0.0,height,height/elevation.shape[0],'d'):
      yCoords.InsertNextValue(y)
    zCoords = vtk.vtkDoubleArray()
    zCoords.InsertNextValue(0.0)

    grid.SetXCoordinates(xCoords)
    grid.SetYCoordinates(yCoords)
    grid.SetZCoordinates(zCoords)

    vector = vtk.vtkFloatArray()
    vector.SetName("elevation");
    vector.SetNumberOfComponents(3);
    vector.SetNumberOfTuples(numpy.prod(elevation.shape));
    for i,z in enumerate(numpy.ravel(elevation)):
      vector.SetTuple3(i,0,0,z)

    grid.GetPointData().AddArray(vector)

    writer = vtk.vtkXMLRectilinearGridWriter()
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    writer.SetFileName(os.path.splitext(file)[0]+'.vtr')
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(grid)
    else:
        writer.SetInputData(grid)
    writer.Write()
