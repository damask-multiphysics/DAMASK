#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,sys,re
import argparse
import damask
import vtk, numpy as np

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName, damask.version])

parser = argparse.ArgumentParser(description='Convert from Marc input file format (.dat) to VTK format (.vtu)', version = scriptID)
parser.add_argument('filename', type=str, help='file to convert')
parser.add_argument('-t', '--table',  type=str, help='ASCIItable file containing nodal data to subdivide and interpolate')

args = parser.parse_args()

with open(args.filename, 'r') as marcfile:
  marctext = marcfile.read();
  
# Load table (if any)
if args.table is not None:
  try:
    table = damask.ASCIItable(
      name=args.table,
      outname='subdivided_{}'.format(args.table),
      buffered=True
    )
    
    table.head_read()
    table.data_readArray()
    
    # Python list is faster for appending
    nodal_data = list(table.data)
  except: args.table = None

# Extract connectivity chunk from file...
connectivity_text = re.findall(r'connectivity[\n\r]+(.*?)[\n\r]+[a-zA-Z]', marctext, flags=(re.MULTILINE | re.DOTALL))[0]
connectivity_lines = re.split(r'[\n\r]+', connectivity_text, flags=(re.MULTILINE | re.DOTALL))
connectivity_header = connectivity_lines[0]
connectivity_lines = connectivity_lines[1:]

# Construct element map
elements = dict(map(lambda line: 
  (
    int(line[0:10]), # index
    { 
      'type':  int(line[10:20]), 
      'verts': list(map(int, re.split(r' +', line[20:].strip())))
    }
  ), connectivity_lines))

# Extract coordinate chunk from file
coordinates_text = re.findall(r'coordinates[\n\r]+(.*?)[\n\r]+[a-zA-Z]', marctext, flags=(re.MULTILINE | re.DOTALL))[0]
coordinates_lines = re.split(r'[\n\r]+', coordinates_text, flags=(re.MULTILINE | re.DOTALL))
coordinates_header = coordinates_lines[0]
coordinates_lines = coordinates_lines[1:]

# marc input file does not use "e" in scientific notation, this adds it and converts
fl_format = lambda string: float(re.sub(r'(\d)([\+\-])', r'\1e\2', string))
# Construct coordinate map
coordinates = dict(map(lambda line:
  (
    int(line[0:10]),
    np.array([
      fl_format(line[10:30]), 
      fl_format(line[30:50]), 
      fl_format(line[50:70])
    ])
  ), coordinates_lines))

# Subdivide volumes
grid = vtk.vtkUnstructuredGrid()
vertex_count = len(coordinates)
edge_to_vert = dict() # when edges are subdivided, a new vertex in the middle is produced and placed in here
ordered_pair = lambda a, b: (a, b) if a < b else (b, a) # edges are bidirectional

def subdivide_edge(vert1, vert2):
  edge = ordered_pair(vert1, vert2)
  
  if edge in edge_to_vert:
    return edge_to_vert[edge]
  
  # Vertex does not exist, create it
  newvert = len(coordinates) + 1
  coordinates[newvert] = 0.5 * (coordinates[vert1] + coordinates[vert2]) # Average
  edge_to_vert[edge] = newvert;
  
  # Interpolate nodal data
  if args.table is not None:
    nodal_data.append(0.5 * (nodal_data[vert1 - 1] + nodal_data[vert2 - 1]))
  return newvert;

for el_id in range(1, len(elements) + 1): # Marc starts counting at 1
  el = elements[el_id]
  if el['type'] == 7:
    # Hexahedron, subdivided
    
    # There may be a better way to iterate over these, but this is consistent
    # with the ordering scheme provided at https://damask.mpie.de/pub/Documentation/ElementType
    
    subverts = np.zeros((3,3,3), dtype=int)
    # Get corners
    subverts[0, 0, 0] = el['verts'][0]
    subverts[2, 0, 0] = el['verts'][1]
    subverts[2, 2, 0] = el['verts'][2]
    subverts[0, 2, 0] = el['verts'][3]
    subverts[0, 0, 2] = el['verts'][4]
    subverts[2, 0, 2] = el['verts'][5]
    subverts[2, 2, 2] = el['verts'][6]
    subverts[0, 2, 2] = el['verts'][7]
    
    # lower edges
    subverts[1, 0, 0] = subdivide_edge(subverts[0, 0, 0], subverts[2, 0, 0])
    subverts[2, 1, 0] = subdivide_edge(subverts[2, 0, 0], subverts[2, 2, 0])
    subverts[1, 2, 0] = subdivide_edge(subverts[2, 2, 0], subverts[0, 2, 0])
    subverts[0, 1, 0] = subdivide_edge(subverts[0, 2, 0], subverts[0, 0, 0])
    
    # middle edges
    subverts[0, 0, 1] = subdivide_edge(subverts[0, 0, 0], subverts[0, 0, 2])
    subverts[2, 0, 1] = subdivide_edge(subverts[2, 0, 0], subverts[2, 0, 2])
    subverts[2, 2, 1] = subdivide_edge(subverts[2, 2, 0], subverts[2, 2, 2])
    subverts[0, 2, 1] = subdivide_edge(subverts[0, 2, 0], subverts[0, 2, 2])
    
    # top edges
    subverts[1, 0, 2] = subdivide_edge(subverts[0, 0, 2], subverts[2, 0, 2])
    subverts[2, 1, 2] = subdivide_edge(subverts[2, 0, 2], subverts[2, 2, 2])
    subverts[1, 2, 2] = subdivide_edge(subverts[2, 2, 2], subverts[0, 2, 2])
    subverts[0, 1, 2] = subdivide_edge(subverts[0, 2, 2], subverts[0, 0, 2])
    
    # then faces...  The edge_to_vert addition is due to there being two ways
    # to calculate a face vertex, depending on which opposite vertices are used to subdivide.
    # This way, we avoid creating duplicate vertices.
    subverts[1, 1, 0] = subdivide_edge(subverts[1, 0, 0], subverts[1, 2, 0])
    edge_to_vert[ordered_pair(subverts[0, 1, 0], subverts[2, 1, 0])] = subverts[1, 1, 0]
    
    subverts[1, 0, 1] = subdivide_edge(subverts[1, 0, 0], subverts[1, 0, 2])
    edge_to_vert[ordered_pair(subverts[0, 0, 1], subverts[2, 0, 1])] = subverts[1, 0, 1]
    
    subverts[2, 1, 1] = subdivide_edge(subverts[2, 1, 0], subverts[2, 1, 2])
    edge_to_vert[ordered_pair(subverts[2, 0, 1], subverts[2, 2, 1])] = subverts[2, 1, 1]
    
    subverts[1, 2, 1] = subdivide_edge(subverts[1, 2, 0], subverts[1, 2, 2])
    edge_to_vert[ordered_pair(subverts[0, 2, 1], subverts[2, 2, 1])] = subverts[1, 2, 1]
    
    subverts[0, 1, 1] = subdivide_edge(subverts[0, 1, 0], subverts[0, 1, 2])
    edge_to_vert[ordered_pair(subverts[0, 0, 1], subverts[0, 2, 1])] = subverts[0, 1, 1]
    
    subverts[1, 1, 2] = subdivide_edge(subverts[1, 0, 2], subverts[1, 2, 2])
    edge_to_vert[ordered_pair(subverts[0, 1, 2], subverts[2, 1, 2])] = subverts[1, 1, 2]
    
    # and finally the center.  There are three ways to calculate, but elements should
    # not intersect, so the edge_to_vert part isn't needed here.
    subverts[1, 1, 1] = subdivide_edge(subverts[1, 1, 0], subverts[1, 1, 2])
    
    
    # Now make the hexahedron subelements
    # order in which vtk expects vertices for a hexahedron
    order = np.array([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)])
    for z in range(2):
      for y in range(2):
        for x in range(2):
          hex_ = vtk.vtkHexahedron()
          for vert_id in range(8):
            coord = order[vert_id] + (x, y, z)
            # minus one, since vtk starts at zero but marc starts at one
            hex_.GetPointIds().SetId(vert_id, subverts[coord[0], coord[1], coord[2]] - 1) 
          grid.InsertNextCell(hex_.GetCellType(), hex_.GetPointIds())
    
  else:
    damask.util.croak('Unsupported Marc element type: {} (skipping)'.format(el['type']))

# Load all points
points = vtk.vtkPoints()
for point in range(1, len(coordinates) + 1): # marc indices start at 1
  points.InsertNextPoint(coordinates[point].tolist())

grid.SetPoints(points)

# grid now contains the elements from the given marc file
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(re.sub(r'\..+', ".vtu", args.filename))  # *.vtk extension does not work in paraview

if vtk.VTK_MAJOR_VERSION <= 5: writer.SetInput(grid)
else:                          writer.SetInputData(grid)
writer.Write()

if args.table is not None:  
  table.info_append([
    scriptID + ' ' + ' '.join(sys.argv[1:]),
  ])
  table.head_write()
  table.output_flush()
  
  table.data = np.array(nodal_data)
  
  table.data_writeArray()

table.close()
