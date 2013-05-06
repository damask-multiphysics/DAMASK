#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
 
from optparse import OptionParser, Option
from vtk import *
import os,numpy,string,math
 

# -----------------------------
def getHeader(filename,sizeFastIndex,sizeSlowIndex,stepsize):
# -----------------------------
# returns header for ang file
# 
  
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
    '# XSTEP: ' + str(stepsize), \
    '# YSTEP: ' + str(stepsize), \
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
# -----------------------------
# returns positive angle in radians
# gets angle in degrees

  angle = math.radians(float(angle))
  while angle < 0.0:
    angle += 2.0 * math.pi

  return angle


# -----------------------------
def getDataLine(angles,position,validData=True):
# -----------------------------
# returns string of one line in ang file
# convention in ang file: y coordinate comes first and is fastest index
  
  info = {True:  (9999.9, 1.0, 0,99999,0.0),
          False: (  -1.0,-1.0,-1,   -1,1.0)}
  return '%9.5f %9.5f %9.5f %12.5f %12.5f %6.1f %6.3f %2i %6i %6.3f \n'%(tuple(map(positiveRadians,angles))+tuple(position[1::-1])+info[validData])



# --------------------------------------------------------------------
# MAIN FUNCTION STARTS HERE
# --------------------------------------------------------------------

parser = OptionParser(usage='%prog options [file[s]]', description = """
Builds a ang files from a vtk file.

""" + string.replace('$Id$','\n','\\n')
)


parser.add_option('--disp','--displacement',dest='dispLabel', type='string', \
                                            help='label of displacements [%default]')
parser.add_option('--euler',                dest='eulerLabel', type='string', nargs=3, \
                                            help='labels of euler angles [%default]')
parser.add_option('-n','--normal',          dest='normal', type='float', nargs=3, \
                                            help='normal of slices in direction of increasing slice numbers [%default]')
parser.add_option('-u','--up',              dest='up', type='float', nargs=3,
                                            help='up direction of slices [%default]')
parser.add_option('-i','--slices',          dest='Nslices', type='int', \
                                            help='number of slices [%default]')
parser.add_option('-d','--distance',        dest='distance', type='float', \
                                            help='slice distance [%default]')
parser.add_option('-s','--size',            dest='size', type='float', nargs=3, \
                                            help='physical size of ang file [%default]')
parser.add_option('-r','--resolution',      dest='resolution', type='float',
                                            help='scaling factor for resolution [%default]')
parser.add_option('--hex','--hexagonal',    dest='hexagonal', action='store_true',
                                            help='use in plane hexagonal grid [%default]')
parser.add_option('--ds','--dispscaling',   dest='dispScaling', type='float', \
                                            help='scaling of displacements [%default]')
parser.add_option('--interpolation',        dest='interpolation', type='int', \
                                            help='number of points for linear interpolation [%default]')
parser.add_option('--verbose',              dest='verbose', action='store_true',
                                            help='verbose mode [%default]')

parser.set_defaults(dispLabel = 'displacement')
parser.set_defaults(eulerLabel = ['1_eulerangles','2_eulerangles','3_eulerangles'])
parser.set_defaults(hexagonal = False)
parser.set_defaults(normal = [0.0,0.0,1.0])
parser.set_defaults(up = [0.0,1.0,0.0])
parser.set_defaults(Nslices = 1)
parser.set_defaults(distance = 0.0)
parser.set_defaults(size = [1.0,1.0,0.0])
parser.set_defaults(resolution = 1.0)
parser.set_defaults(dispScaling = 1.0)
parser.set_defaults(verbose = False)
parser.set_defaults(interpolation = 1)
(options,filenames) = parser.parse_args()


#--- SANITY CHECKS

# check for valid filenames

for filename in filenames:
  if not os.path.exists(filename):
    parser.error('file "%s" does not exist'%filename)
  if not os.path.splitext(filename)[1] == '.vtk':
    parser.error('"%s": need vtk file'%filename)


# check for othogonality of normal and up vector

if numpy.dot(numpy.array(options.normal),numpy.array(options.up)) > 1e-3:
  parser.error('normal vector and up vector have to be orthogonal')


# check for options that are not yet implemented

if options.interpolation > 1:
  parser.error('interpolation not yet supported')
if options.hexagonal:
  parser.error('hexagonal grid not yet supported')



#--- ITERATE OVER FILES AND PROCESS THEM

for filename in filenames:
  
  # Read the source file
  
  reader = vtkUnstructuredGridReader()
  reader.SetFileName(filename)
  reader.ReadAllScalarsOn()
  reader.ReadAllVectorsOn()
  reader.Update() 
  undeformedMesh = reader.GetOutput()
  
  
  # Get euler angles from cell data
  
  angles = {}
  Nscalars = reader.GetNumberOfScalarsInFile()
  for i in range(Nscalars):
    scalarName = reader.GetScalarsNameInFile(i)
    if scalarName in options.eulerLabel:
      angles[scalarName] = undeformedMesh.GetCellData().GetScalars(scalarName)
  if len(angles) < 3:    # found data for all three euler angles?
    for label in options.eulerLabel:
      if not label in angles.keys():
        parser.error('Could not find scalar data with name %s'%label)
  
  
  # Get deformed mesh
  
  warpVector = vtkWarpVector()
  warpVector.SetScaleFactor(options.dispScaling)
  warpVector.SetInput(undeformedMesh)
  warpVector.Update()
  deformedMesh = warpVector.GetOutput()       # todo: not clear how to choose other vector data than the first entry
  box = deformedMesh.GetBounds()              # bounding box in mesh system
  if options.verbose:
    print ''
    print 'MESH SYSTEM'
    print '  bounding box'
    print '    x ',[box[0],box[1]]
    print '    y ',[box[2],box[3]]
    print '    z ',[box[4],box[5]]


  # Get cell centers of deformed mesh (position of ips)
  
  cellCenter = vtkCellCenters()
  cellCenter.SetVertexCells(0)   # do not generate vertex cells, just points
  cellCenter.SetInput(deformedMesh)
  cellCenter.Update()
  meshIPs = cellCenter.GetOutput()


  # Get outer surface of deformed mesh
  
  surfaceFilter = vtkDataSetSurfaceFilter()
  surfaceFilter.SetInput(deformedMesh)
  surfaceFilter.Update()
  surface = surfaceFilter.GetOutput()
  
  
  # Get coordinate system for ang files
  # z-vector is normal to slices
  # x-vector corresponds to the up-direction
  # "R" rotates coordinates from the mesh system into the TSL system

  z = numpy.array(options.normal,dtype='float')
  z = z / numpy.linalg.norm(z)
  x = numpy.array(options.up,dtype='float')
  x = x / numpy.linalg.norm(x)
  y = numpy.cross(z,x)
  R = numpy.array([x,y,z])


  # Get bounding box in rotated system (x,y,z)

  rotatedbox = [[numpy.inf,-numpy.inf] for i in range(3)]     # bounding box in rotated TSL system
  for n in range(8):                                          # loop over eight vertices of mesh bounding box 
    vert = numpy.array([box[0+(n/1)%2], 
                        box[2+(n/2)%2], 
                        box[4+(n/4)%2]])                      # vertex in mesh system
    rotatedvert = numpy.dot(R,vert)                           # vertex in rotated system
    for i in range(3):
      rotatedbox[i][0] = min(rotatedbox[i][0],rotatedvert[i])
      rotatedbox[i][1] = max(rotatedbox[i][1],rotatedvert[i])


  # Correct bounding box so that a multiplicity of the resolution fits into it
  # and get number of points and extent in each (rotated) axis direction

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
  if options.verbose:
    print ''
    print 'ROTATED SYSTEM'
    print '  axis (x: up direction, z: slice normal)'
    print '    x ',list(x)
    print '    y ',list(y)
    print '    z ',list(z)
    print '  bounding box'
    print '    x ',rotatedbox[0]
    print '    y ',rotatedbox[1]
    print '    z ',rotatedbox[2]
    print '  number of points per slice'
    print '    x ',Npoints[0]
    print '    y ',Npoints[1]
    print '  number of slices'
    print '    z ',Npoints[2]


  # Generate new regular point grid for ang files
  # Use "polydata" object with points as single vertices
  # beware of TSL convention: y direction is fastest index

  points = vtkPoints()
  for k in xrange(Npoints[2]):
    for j in xrange(Npoints[0]):  
      for i in xrange(Npoints[1]):   # y is fastest index
        rotatedpoint = numpy.array([rotatedbox[0][0] + (float(j) + 0.5) * options.resolution,
                                    rotatedbox[1][0] + (float(i) + 0.5) * options.resolution,
                                    rotatedbox[2][0] + (float(k) + 0.5) * options.distance ])  # point in rotated system
        point = numpy.dot(R.T,rotatedpoint)                                                    # point in mesh system
        points.InsertNextPoint(list(point))

  vertices = vtkCellArray()
  for i in xrange(Npoints[0]*Npoints[1]*Npoints[2]):
    vertex = vtkVertex()
    vertex.GetPointIds().SetId(0,i)  # each vertex consists of exactly one (index 0) point with ID "i"
    vertices.InsertNextCell(vertex)

  pointgrid = vtk.vtkPolyData()
  pointgrid.SetPoints(points)
  pointgrid.SetVerts(vertices)
  pointgrid.Update()


  # Find out which points reside inside mesh geometry

  enclosedPoints = vtkSelectEnclosedPoints()
  enclosedPoints.SetSurface(surface)
  enclosedPoints.SetInput(pointgrid)
  enclosedPoints.Update()


  # Build kdtree from mesh IPs and match mesh IPs to point grid
  # could also be done with nearest neighbor search from damask.core, would possibly be faster ?

  kdTree = vtkKdTree()
  kdTree.BuildLocatorFromPoints(meshIPs.GetPoints())
  gridToMesh = []
  ids = vtkIdList()
  for i in range(pointgrid.GetNumberOfPoints()):
    gridToMesh.append([])
    if enclosedPoints.IsInside(i):
      kdTree.FindClosestNPoints(options.interpolation,pointgrid.GetPoint(i),ids) # here one could use faster(?) "FindClosestPoint" if only first nearest neighbor required
      for j in range(ids.GetNumberOfIds()): 
        gridToMesh[-1].extend([ids.GetId(j)])
  


  # ITERATE OVER SLICES AND CREATE ANG FILE
  
  NpointsPerSlice = Npoints[0] * Npoints[1]
  for sliceN in range(Npoints[2]):
    
    # Open file and write header
    
    angfilename = eval('"'+eval("'%%s_slice%%0%ii.ang'%(math.log10(Npoints[2])+1)")+'"%(os.path.splitext(filename)[0],sliceN)')
    with open(angfilename,'w') as angfile:
      angfile.write(getHeader(filename,Npoints[1],Npoints[0],options.resolution))
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

        angfile.write(getDataLine(interpolatedPhi,pointgrid.GetPoint(i),enclosedPoints.IsInside(i)))
      
