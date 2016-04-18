#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,string,math,sys
import numpy as np
from optparse import OptionParser
import vtk
import damask
 
scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# -----------------------------
def getHeader(filename,sizeFastIndex,sizeSlowIndex,stepsize):
  """returns header for ang file step size in micrometer"""
  return '\n'.join([ \
    '# TEM_PIXperUM          1.000000', \
    '# x-star                1.000000', \
    '# y-star                1.000000', \
    '# z-star                1.000000', \
    '# WorkingDistance       18.000000', \
    '#', \
    '# Phase                 1', \
    '# MaterialName          XX', \
    '# Formula               XX', \
    '# Info', \
    '# Symmetry              43', \
    '# LatticeConstants      2.870 2.870 2.870  90.000  90.000  90.000', \
    '# NumberFamilies        1', \
    '# hklFamilies           1  1  0 1 0.000000 1', \
    '# Categories            0 0 0 0 0 ', \
    '#', \
    '# GRID: SqrGrid', \
    '# XSTEP: ' + str(stepsize*1e6), \
    '# YSTEP: ' + str(stepsize*1e6), \
    '# NCOLS_ODD: ' + str(sizeFastIndex), \
    '# NCOLS_EVEN: ' + str(sizeFastIndex), \
    '# NROWS: ' + str(sizeSlowIndex), \
    '#', \
    '# OPERATOR: ' + string.replace('$Id$','\n','\\n'), \
    '#', \
    '# SAMPLEID: %s'%filename, \
    '#', \
    '# SCANID: ', \
    '#', \
    ]) + '\n'


# -----------------------------
def positiveRadians(angle):
  """returns positive angle in radians from angle in degrees"""
  angle = math.radians(float(angle))
  while angle < 0.0:
    angle += 2.0 * math.pi

  return angle


# -----------------------------
def getDataLine(angles,x,y,validData=True):
  """
  returns string of one line in ang file

  convention in ang file: y coordinate comes first and is fastest index
  positions in micrometer
  """
  info = {True:  (9999.9, 1.0, 0,99999,0.0),
          False: (  -1.0,-1.0,-1,   -1,1.0)}
  return '%9.5f %9.5f %9.5f %12.5f %12.5f %6.1f %6.3f %2i %6i %6.3f \n'\
                %(tuple(map(positiveRadians,angles))+(y*1e6,x*1e6)+info[validData])



# --------------------------------------------------------------------
# MAIN FUNCTION STARTS HERE
# --------------------------------------------------------------------

parser = OptionParser(usage='%prog options [file[s]]', description = """
Builds a ang files from a vtk file.

""", version = scriptID)


parser.add_option('--disp','--displacement',dest='dispLabel', \
                                            help='label of displacements [%default]')
parser.add_option('--euler',                dest='eulerLabel', nargs=3, \
                                            help='labels of euler angles [%default]')
parser.add_option('-n','--normal',          dest='normal', type='float', nargs=3, \
                                            help='normal of slices in direction of increasing slice numbers [%default]')
parser.add_option('-u','--up',              dest='up', type='float', nargs=3,
                                            help='up direction of slices [%default]')
parser.add_option('-i','--slices',          dest='Nslices', type='int', \
                                            help='number of slices [%default]')
parser.add_option('-d','--distance',        dest='distance', type='float', \
                                            help='slice distance [%default]')
parser.add_option('-s','--scale',           dest='scale', type='float', \
                                            help='scale length from vtk file [%default]')
parser.add_option('-r','--resolution',      dest='resolution', type='float',
                                            help='scaling factor for resolution [%default]')
parser.add_option('--hex','--hexagonal',    dest='hexagonal', action='store_true',
                                            help='use in plane hexagonal grid [%default]')
parser.add_option('--interpolation',        dest='interpolation', type='int', \
                                            help='number of points for linear interpolation [%default]')
parser.add_option('--verbose',              dest='verbose', action='store_true',
                                            help='verbose mode [%default]')
parser.add_option('--visualize',            dest='visualize', action='store_true',
                                            help='visualize geometry [%default]')

parser.set_defaults(dispLabel = 'displacement')
parser.set_defaults(eulerLabel = ['1_1_eulerangles','1_2_eulerangles','1_3_eulerangles'])
parser.set_defaults(hexagonal = False)
parser.set_defaults(normal = [0.0,0.0,-1.0])
parser.set_defaults(up = [0.0,1.0,0.0])
parser.set_defaults(Nslices = 1)
parser.set_defaults(distance = 0.0)
parser.set_defaults(scale = 1.0)
parser.set_defaults(resolution = 1.0)
parser.set_defaults(dispScaling = 1.0)
parser.set_defaults(interpolation = 1)
parser.set_defaults(verbose = False)
parser.set_defaults(visualize = False)
(options,filenames) = parser.parse_args()


#--- SANITY CHECKS

# check for valid filenames

for filename in filenames:
  if not os.path.exists(filename):
    parser.error('file "%s" does not exist'%filename)
  if not os.path.splitext(filename)[1] == '.vtk':
    parser.error('"%s": need vtk file'%filename)


# check for othogonality of normal and up vector

if np.dot(np.array(options.normal),np.array(options.up)) > 1e-3:
  parser.error('normal vector and up vector have to be orthogonal')


# check for options that are not yet implemented

if options.interpolation > 1:
  parser.error('interpolation not yet supported')
if options.hexagonal:
  parser.error('hexagonal grid not yet supported')



#--- ITERATE OVER FILES AND PROCESS THEM

for filename in filenames:
  
  if options.verbose: sys.stdout.write("\nREADING VTK FILE\n")
# Read the source file
  reader = vtk.vtkUnstructuredGridReader()
  reader.SetFileName(filename)
  reader.ReadAllScalarsOn()
  reader.ReadAllVectorsOn()
  reader.Update() 
  undeformedMesh = reader.GetOutput()

  
# Get euler angles from cell data
  
  if options.verbose: sys.stdout.write("\nGETTING EULER ANGLES\n")
  angles = {}
  for i in range(reader.GetNumberOfScalarsInFile()):
    scalarName = reader.GetScalarsNameInFile(i)
    if scalarName in options.eulerLabel:
      angles[scalarName] = undeformedMesh.GetCellData().GetScalars(scalarName)
      if options.verbose: sys.stdout.write("  found scalar with name %s\n"%scalarName)
  if len(angles) < 3:    # found data for all three euler angles?
    for label in options.eulerLabel:
      if label not in angles.keys():
        parser.error('Could not find scalar data with name %s'%label)
  
  
# Get deformed mesh
  
  if options.verbose: sys.stdout.write("\nDEFORM MESH\n")
  warpVector = vtk.vtkWarpVector()
  undeformedMesh.GetPointData().SetActiveVectors(options.dispLabel)
  warpVector.SetInput(undeformedMesh)
  warpVector.Update()
  deformedMesh = warpVector.GetOutput()
  box = deformedMesh.GetBounds()              # bounding box in mesh system
  if options.verbose:
    sys.stdout.write("  bounding box in lab system\n")
    sys.stdout.write("    x (% .8f % .8f)\n"%(box[0],box[1]))
    sys.stdout.write("    y (% .8f % .8f)\n"%(box[2],box[3]))
    sys.stdout.write("    z (% .8f % .8f)\n"%(box[4],box[5]))


# Get cell centers of deformed mesh (position of ips)
  
  if options.verbose: sys.stdout.write("\nGETTING CELL CENTERS OF DEFORMED MESH\n")
  cellCenter = vtk.vtkCellCenters()
  cellCenter.SetVertexCells(0)   # do not generate vertex cells, just points
  cellCenter.SetInput(deformedMesh)
  cellCenter.Update()
  meshIPs = cellCenter.GetOutput()


# Get outer surface of deformed mesh
  
  if options.verbose: sys.stdout.write("\nGETTING OUTER SURFACE OF DEFORMED MESH\n")
  surfaceFilter = vtk.vtkDataSetSurfaceFilter()
  surfaceFilter.SetInput(deformedMesh)
  surfaceFilter.Update()
  surface = surfaceFilter.GetOutput()
  
  
# Get coordinate system for ang files
# z-vector is normal to slices
# x-vector corresponds to the up-direction
# "R" rotates coordinates from the mesh system into the TSL system

  if options.verbose: sys.stdout.write("\nGETTING COORDINATE SYSTEM FOR ANG FILES\n")
  z = np.array(options.normal,dtype='float')
  z = z / np.linalg.norm(z)
  x = np.array(options.up,dtype='float')
  x = x / np.linalg.norm(x)
  y = np.cross(z,x)
  R = np.array([x,y,z])
  if options.verbose:
    sys.stdout.write("  axis (x: up direction, z: slice normal)\n")
    sys.stdout.write("    x (% .8f % .8f % .8f)\n"%tuple(x))
    sys.stdout.write("    y (% .8f % .8f % .8f)\n"%tuple(y))
    sys.stdout.write("    z (% .8f % .8f % .8f)\n"%tuple(z))


# Get bounding box in rotated system (x,y,z)

  if options.verbose: sys.stdout.write("\nGETTING BOUNDING BOX IN ROTATED SYSTEM\n")
  rotatedbox = [[np.inf,-np.inf] for i in range(3)]     # bounding box in rotated TSL system
  for n in range(8):                                          # loop over eight vertices of mesh bounding box 
    vert = np.array([box[0+(n/1)%2], 
                        box[2+(n/2)%2], 
                        box[4+(n/4)%2]])                      # vertex in mesh system
    rotatedvert = np.dot(R,vert)                           # vertex in rotated system
    for i in range(3):
      rotatedbox[i][0] = min(rotatedbox[i][0],rotatedvert[i])
      rotatedbox[i][1] = max(rotatedbox[i][1],rotatedvert[i])
  if options.verbose:
    sys.stdout.write("  bounding box in rotated system\n")
    sys.stdout.write("    x (% .8f % .8f)\n"%tuple(rotatedbox[0]))
    sys.stdout.write("    y (% .8f % .8f)\n"%tuple(rotatedbox[1]))
    sys.stdout.write("    z (% .8f % .8f)\n"%tuple(rotatedbox[2]))


# Correct bounding box so that a multiplicity of the resolution fits into it
# and get number of points and extent in each (rotated) axis direction

  if options.verbose: sys.stdout.write("\nCORRECTING EXTENT OF BOUNDING BOX IN ROTATED SYSTEM\n")
  correction = []
  Npoints = []
  extent = [rotatedbox[i][1] - rotatedbox[i][0] for i in range(3)]
  for i in range(2):
    Npoints.extend([int(math.ceil(extent[i] / options.resolution))])
    correction.extend([float(Npoints[i]) * options.resolution - extent[i]])
  if options.distance > 0.0: 
    Npoints.extend([int(math.ceil(extent[2] / options.distance))])
    correction.extend([float(Npoints[2]) * options.distance - extent[2]])
  else:
    Npoints.extend([options.Nslices])
    correction.extend([0.0])
    options.distance = extent[2] / float(options.Nslices)
  for i in range(3):
    rotatedbox[i][0] = rotatedbox[i][0] - 0.5 * correction[i]
    rotatedbox[i][1] = rotatedbox[i][1] + 0.5 * correction[i]
    extent[i] = rotatedbox[i][1] - rotatedbox[i][0]
  NpointsPerSlice = Npoints[0] * Npoints[1]
  totalNpoints = NpointsPerSlice * Npoints[2]
  if options.verbose:
    sys.stdout.write("  corrected bounding box in rotated system\n")
    sys.stdout.write("    x (% .8f % .8f)\n"%tuple(rotatedbox[0]))
    sys.stdout.write("    y (% .8f % .8f)\n"%tuple(rotatedbox[1]))
    sys.stdout.write("    z (% .8f % .8f)\n"%tuple(rotatedbox[2]))


# Generate new regular point grid for ang files
# Use "polydata" object with points as single vertices
# beware of TSL convention: y direction is fastest index

  if options.verbose: sys.stdout.write("\nGENERATING POINTS FOR POINT GRID")
  points = vtk.vtkPoints()
  for k in xrange(Npoints[2]):
    for j in xrange(Npoints[0]):  
      for i in xrange(Npoints[1]):   # y is fastest index
        rotatedpoint = np.array([rotatedbox[0][0] + (float(j) + 0.5) * options.resolution,
                                    rotatedbox[1][0] + (float(i) + 0.5) * options.resolution,
                                    rotatedbox[2][0] + (float(k) + 0.5) * options.distance ])  # point in rotated system
        point = np.dot(R.T,rotatedpoint)                                                    # point in mesh system
        points.InsertNextPoint(list(point))
        if options.verbose: 
          sys.stdout.write("\rGENERATING POINTS FOR POINT GRID %d%%" %(100*(Npoints[1]*(k*Npoints[0]+j)+i+1)/totalNpoints))
          sys.stdout.flush()
  if options.verbose: 
    sys.stdout.write("\n  number of slices: %i\n"%Npoints[2])
    sys.stdout.write("  slice spacing: %.8f\n"%options.distance)
    if Npoints[2] > 1: 
      sys.stdout.write("  number of points per slice: %i = %i rows * %i points in row\n"%(NpointsPerSlice,Npoints[0],Npoints[1]))
    sys.stdout.write("  grid resolution: %.8f\n"%options.resolution)

  if options.verbose: sys.stdout.write("\nGENERATING VERTICES FOR POINT GRID")
  vertices = vtk.vtkCellArray()
  for i in xrange(totalNpoints):
    vertex = vtk.vtkVertex()
    vertex.GetPointIds().SetId(0,i)  # each vertex consists of exactly one (index 0) point with ID "i"
    vertices.InsertNextCell(vertex)
    if options.verbose: 
      sys.stdout.write("\rGENERATING VERTICES FOR POINT GRID %d%%" %(100*(i+1)/totalNpoints))
      sys.stdout.flush()

  if options.verbose: sys.stdout.write("\n\nGENERATING POINT GRID\n")
  pointgrid = vtk.vtkPolyData()
  pointgrid.SetPoints(points)
  pointgrid.SetVerts(vertices)
  pointgrid.Update()


# Find out which points reside inside mesh geometry

  if options.verbose: sys.stdout.write("\nIDENTIFYING POINTS INSIDE MESH GEOMETRY\n")
  enclosedPoints = vtk.vtkSelectEnclosedPoints()
  enclosedPoints.SetSurface(surface)
  enclosedPoints.SetInput(pointgrid)
  enclosedPoints.Update()


# Build kdtree from mesh IPs and match mesh IPs to point grid

  if options.verbose: sys.stdout.write("\nBUILDING MAPPING OF GRID POINTS")
  kdTree = vtk.vtkKdTree()
  kdTree.BuildLocatorFromPoints(meshIPs.GetPoints())
  gridToMesh = []
  ids = vtk.vtkIdList()
  NenclosedPoints = 0
  for i in range(pointgrid.GetNumberOfPoints()):
    gridToMesh.append([])
    if enclosedPoints.IsInside(i):
      NenclosedPoints += 1
# here one could use faster(?) "FindClosestPoint" if only first nearest neighbor required
      kdTree.FindClosestNPoints(options.interpolation,pointgrid.GetPoint(i),ids) 
      for j in range(ids.GetNumberOfIds()): 
        gridToMesh[-1].extend([ids.GetId(j)])
    if options.verbose: 
      sys.stdout.write("\rBUILDING MAPPING OF GRID POINTS %d%%" %(100*(i+1)/totalNpoints))
      sys.stdout.flush()
  if options.verbose:
    sys.stdout.write("\n  Number of points inside mesh geometry %i\n"%NenclosedPoints)
    sys.stdout.write("  Number of points outside mesh geometry %i\n"%(totalNpoints - NenclosedPoints))
  


# ITERATE OVER SLICES AND CREATE ANG FILE
  
  if options.verbose: 
    sys.stdout.write("\nWRITING OUT ANG FILES\n")
    sys.stdout.write("  scaling all length with %f\n"%options.scale)
  x0,y0,z0 = np.dot(R,pointgrid.GetPoint(0))                    # first point on slice defines origin
  for sliceN in range(Npoints[2]):
    
    # Open file and write header
    
    angfilename = eval('"'+eval("'%%s_slice%%0%ii.ang'%(math.log10(Npoints[2])+1)")+'"%(os.path.splitext(filename)[0],sliceN+1)')
    with open(angfilename,'w') as angfile:
      if options.verbose: sys.stdout.write("  %s\n"%angfilename)
      angfile.write(getHeader(filename,Npoints[1],Npoints[0],options.resolution*options.scale))
      for i in xrange(sliceN*NpointsPerSlice,(sliceN+1)*NpointsPerSlice):   # Iterate over points on slice
        

        # Get euler angles of closest IDs

        if enclosedPoints.IsInside(i):
          phi = []
          for j in range(len(gridToMesh[i])):
            IP = gridToMesh[i][j]
            phi.append([])
            for k in range(3):
              phi[-1].extend([angles[options.eulerLabel[k]].GetValue(IP)])
        else:
          phi = [[720,720,720]] # fake angles
        
        
        # Interpolate Euler angle
        # NOT YET IMPLEMENTED, simply take the nearest neighbors values

        interpolatedPhi = phi[0]
        
        
        # write data to ang file

        x,y,z = np.dot(R,pointgrid.GetPoint(i))                  # point in rotated TSL system
        x -= x0                                                     # first point on slice defines origin 
        y -= y0                                                     # first point on slice defines origin 
        x *= options.scale
        y *= options.scale
        angfile.write(getDataLine(interpolatedPhi,x,y,enclosedPoints.IsInside(i)))
      

# Visualize slices
  
  if options.visualize:
    meshMapper = vtk.vtkDataSetMapper()
    meshMapper.SetInput(surface)
    meshMapper.ScalarVisibilityOff()          # do not use scalar data for coloring
    meshActor = vtk.vtkActor()
    meshActor.SetMapper(meshMapper)
    meshActor.GetProperty().SetOpacity(0.2)
    meshActor.GetProperty().SetColor(1.0,1.0,0)
    meshActor.GetProperty().BackfaceCullingOn()
    # meshActor.GetProperty().SetEdgeColor(1,1,0.5)
    # meshActor.GetProperty().EdgeVisibilityOn()
    
    boxpoints = vtk.vtkPoints()
    for n in range(8):
      P = [rotatedbox[0][(n/1)%2],
           rotatedbox[1][(n/2)%2], 
           rotatedbox[2][(n/4)%2]]
      boxpoints.InsertNextPoint(list(np.dot(R.T,np.array(P))))
    box = vtk.vtkHexahedron()
    for n,i in enumerate([0,1,3,2,4,5,7,6]):
      box.GetPointIds().SetId(n,i)
    boxgrid = vtk.vtkUnstructuredGrid()
    boxgrid.SetPoints(boxpoints)
    boxgrid.InsertNextCell(box.GetCellType(), box.GetPointIds())
    boxsurfaceFilter = vtk.vtkDataSetSurfaceFilter()
    boxsurfaceFilter.SetInput(boxgrid)
    boxsurfaceFilter.Update()
    boxsurface = boxsurfaceFilter.GetOutput()

    boxMapper = vtk.vtkDataSetMapper()
    boxMapper.SetInput(boxsurface)
    boxActor = vtk.vtkActor()
    boxActor.SetMapper(boxMapper)
    boxActor.GetProperty().SetLineWidth(2.0)
    boxActor.GetProperty().SetRepresentationToWireframe()

    gridMapper = vtk.vtkDataSetMapper()
    gridMapper.SetInput(pointgrid)
    gridActor = vtk.vtkActor()
    gridActor.SetMapper(gridMapper)
    gridActor.GetProperty().SetColor(0,0,0)
    gridActor.GetProperty().SetPointSize(3)


    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.FullScreenOn()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderer.AddActor(meshActor)
    renderer.AddActor(boxActor)
    renderer.AddActor(gridActor)
    renderer.SetBackground(1,1,1)
     
    renderWindow.Render()
    renderWindowInteractor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
    renderWindowInteractor.Start()

