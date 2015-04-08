#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,threading,re,time,string,fnmatch,vtk
import numpy as np
from optparse import OptionParser
from vtk.util import numpy_support
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

def outStdout(cmd,locals):
  if cmd[0:3] == '(!)':
    exec(cmd[3:])
  elif cmd[0:3] == '(?)':
    cmd = eval(cmd[3:])
    print cmd
  else:
    print cmd
  return

def outFile(cmd,locals):
  if cmd[0:3] == '(!)':
    exec(cmd[3:])
  elif cmd[0:3] == '(?)':
    cmd = eval(cmd[3:])
    locals['filepointer'].write(cmd+'\n')
  else:
    locals['filepointer'].write(cmd+'\n')
  return


def output(cmds,locals,dest):
  for cmd in cmds:
    if isinstance(cmd,list):
      output(cmd,locals,dest)
    else:
      {\
      'File': outFile,\
      'Stdout': outStdout,\
      }[dest](str(cmd),locals)
  return


def transliterateToFloat(x):
  try:
    return float(x)
  except:
    return 0.0


def unravel(item):
  if hasattr(item,'__contains__'): return ' '.join(map(unravel,item))
  else: return str(item)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++
def vtk_writeASCII_mesh(mesh,data,res,sep):
# ++++++++++++++++++++++++++++++++++++++++++++++++++++
  """ function writes data array defined on a hexahedral mesh (geometry) """
  info =  {\
             'tensor': {'name':'tensor','len':9},\
             'vector': {'name':'vector','len':3},\
             'scalar': {'name':'scalar','len':1},\
             'double': {'name':'scalar','len':2},\
             'triple': {'name':'scalar','len':3},\
             'quadruple': {'name':'scalar','len':4},\
            }
  N1 = (res[0]+1)*(res[1]+1)*(res[2]+1)
  N  = res[0]*res[1]*res[2]
  
  cmds = [\
          '# vtk DataFile Version 3.1',
          'powered by %s'%scriptID,
          'ASCII',
          'DATASET UNSTRUCTURED_GRID',
          'POINTS %i double'%N1,
          [[['\t'.join(map(str,mesh[:,i,j,k])) for i in range(res[0]+1)] for j in range(res[1]+1)] for k in range(res[2]+1)],
          'CELLS %i %i'%(N,N*9),
          ]

# cells
  for z in range (res[2]):
    for y in range (res[1]):
      for x in range (res[0]):
        base = z*(res[1]+1)*(res[0]+1)+y*(res[0]+1)+x
        cmds.append('8 '+'\t'.join(map(str,[ \
                                            base,
                                            base+1,
                                            base+res[0]+2,
                                            base+res[0]+1,
                                            base+(res[1]+1)*(res[0]+1),
                                            base+(res[1]+1)*(res[0]+1)+1,
                                            base+(res[1]+1)*(res[0]+1)+res[0]+2,
                                            base+(res[1]+1)*(res[0]+1)+res[0]+1,
                                          ])))
  cmds += [\
           'CELL_TYPES %i'%N,
           ['12']*N,
           'CELL_DATA %i'%N,
          ]
  
  for type in data:
    plural = {True:'',False:'S'}[type.lower().endswith('s')]
    for item in data[type]['_order_']:
      cmds += [\
               '%s %s double'%(info[type]['name'].upper()+plural,item),
               {True:'LOOKUP_TABLE default',False:''}[info[type]['name'][:3]=='sca'],
               [[[sep.join(map(unravel,data[type][item][:,j,k]))] for j in range(res[1])] for k in range(res[2])],
              ]

  return cmds
 
# +++++++++++++++++++++++++++++++++++++++++++++++++++
def vtk_writeASCII_points(coordinates,data,res,sep):
# +++++++++++++++++++++++++++++++++++++++++++++++++++
  """ function writes data array defined on a point field """
  N  = res[0]*res[1]*res[2]
  
  cmds = [\
          '# vtk DataFile Version 3.1',
          'powered by %s'%scriptID,
          'ASCII',
          'DATASET UNSTRUCTURED_GRID',
          'POINTS %i double'%N,
          [[['\t'.join(map(str,coordinates[i,j,k])) for i in range(res[0])] for j in range(res[1])] for k in range(res[2])],
          'CELLS %i %i'%(N,N*2),
          ['1\t%i'%i for i in range(N)],
          'CELL_TYPES %i'%N,
          ['1']*N,
          'POINT_DATA %i'%N,
         ]
  
  for type in data:
    plural = {True:'',False:'S'}[type.lower().endswith('s')]
    for item in data[type]:
      cmds += [\
               '%s %s double'%(type.upper()+plural,item),
               {True:'LOOKUP_TABLE default',False:''}[type.lower()[:3]=='sca'],
               [[[sep.join(map(unravel,data[type][item][:,j,k]))] for j in range(res[1])] for k in range(res[2])],
              ]

  return cmds


# ----------------------- MAIN -------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [options] datafile[s]', description = """
Produce VTK file from data field.
Coordinates are taken from (consecutive) x, y, and z columns.

""", version = scriptID)

sepChoices = ['n','t','s']
parser.add_option('-s', '--scalar', dest='scalar', action='extend', metavar = '<string LIST>',
                                    help='list of single scalars to visualize')
parser.add_option(      '--double', dest='double', action='extend', metavar = '<string LIST>',
                                    help='list of two scalars to visualize')
parser.add_option(      '--triple', dest='triple', action='extend', metavar = '<string LIST>',
                                    help='list of three scalars to visualize')
parser.add_option(      '--quadruple', dest='quadruple', action='extend', metavar = '<string LIST>',
                                    help='list of four scalars to visualize')
parser.add_option('-v', '--vector', dest='vector', action='extend', metavar = '<string LIST>', 
                                    help='list of vectors to visualize')
parser.add_option('-t', '--tensor', dest='tensor', action='extend', metavar = '<string LIST>',
                                    help='list of tensors to visualize')
parser.add_option('-d', '--deformation', dest='defgrad', metavar = 'string',
                                    help='heading of deformation gradient columns [%default]')
parser.add_option('--reference',    dest='undeformed', action='store_true',
                                    help='map results to reference (undeformed) configuration [%default]')
parser.add_option('-c','--cell',    dest='cell', action='store_true',
                                    help='data is cell-centered [%default]')
parser.add_option('-p','--vertex',  dest='cell', action='store_false',
                                    help='data is vertex-centered')
parser.add_option('--mesh',         dest='output_mesh', action='store_true',
                                    help='produce VTK mesh file [%default]')
parser.add_option('--nomesh',       dest='output_mesh', action='store_false',
                                    help='omit VTK mesh file')
parser.add_option('--points',       dest='output_points', action='store_true',
                                    help='produce VTK points file [%default]')
parser.add_option('--nopoints',     dest='output_points', action='store_false',
                                    help='omit VTK points file')
parser.add_option('--separator',    dest='separator', type='choice', choices=sepChoices, metavar='string',
                                    help='data separator {%s} [t]'%(' '.join(map(str,sepChoices))))
parser.add_option('--scaling',      dest='scaling', action='extend', metavar = '<float LIST>',
                                    help='scaling of fluctuation')
parser.add_option('-u', '--unitlength', dest='unitlength', type='float', metavar = 'float',
                                    help='set unit length for 2D model [%default]')

parser.set_defaults(defgrad = 'f')
parser.set_defaults(separator = 't')
parser.set_defaults(scalar = [])
parser.set_defaults(double = [])
parser.set_defaults(triple = [])
parser.set_defaults(quadruple = [])
parser.set_defaults(vector = [])
parser.set_defaults(tensor = [])
parser.set_defaults(output_mesh = True)
parser.set_defaults(output_points = False)
parser.set_defaults(scaling = [])
parser.set_defaults(undeformed = False)
parser.set_defaults(unitlength = 0.0)
parser.set_defaults(cell = True)

sep = {'n': '\n', 't': '\t', 's': ' '}

(options, args) = parser.parse_args()

options.scaling += [1.0 for i in xrange(max(0,3-len(options.scaling)))]
options.scaling = map(float, options.scaling)

for filename in args:
  if not os.path.exists(filename):
    continue
  file = open(filename)
  content = file.readlines()
  file.close()
  m = re.search('(\d+)\s*head', content[0].lower())
  if m == None:
    continue
  print filename,'\n'
  sys.stdout.flush()
  
  headrow = int(m.group(1))
  headings = content[headrow].split()
  column = {}
  matches = {}
  maxcol = 0
  locol = -1
  
  for col,head in enumerate(headings):
    if head == {True:'1_ipinitialcoord',False:'1_nodeinitialcoord'}[options.cell]:
      locol = col
      maxcol = max(maxcol,col+3)
      break

  if locol < 0:
    print 'missing coordinates..!'
    continue

  column['tensor'] = {}
  matches['tensor'] = {}
  for label in [options.defgrad] + options.tensor:
    column['tensor'][label] = -1
    for col,head in enumerate(headings):
      if head == label or head == '1_'+label:
        column['tensor'][label] = col
        maxcol = max(maxcol,col+9)
        matches['tensor'][label]  = [label]
        break
      
  if not options.undeformed and column['tensor'][options.defgrad] < 0:
    print 'missing deformation gradient "%s"..!'%options.defgrad
    continue

  column['vector'] = {}
  matches['vector'] = {}
  for label in options.vector:
    column['vector'][label] = -1
    for col,head in enumerate(headings):
      if head == label or head == '1_'+label:
        column['vector'][label] = col
        maxcol = max(maxcol,col+3)
        matches['vector'][label]  = [label]
        break

  for length,what in enumerate(['scalar','double','triple','quadruple']):
    column[what] = {}
    labels = eval("options.%s"%what)
    matches[what] = {}
    for col,head in enumerate(headings):
      for needle in labels:
        if fnmatch.fnmatch(head,needle):
          column[what][head] = col
          maxcol = max(maxcol,col+1+length)
          if needle not in matches[what]:
            matches[what][needle]  = [head]
          else:
            matches[what][needle] += [head]


  values = np.array(sorted([map(transliterateToFloat,line.split()[:maxcol]) for line in content[headrow+1:]],
                              key=lambda x:(x[locol+0],x[locol+1],x[locol+2])),'d')                         # sort with z as fastest and x as slowest index
  values2 = np.array([map(transliterateToFloat,line.split()[:maxcol]) for line in content[headrow+1:]],'d') # sort with x as fastest and z as slowest index

  N = len(values)

  tempGrid = [{},{},{}]
  for j in xrange(3):
    for i in xrange(N):
      tempGrid[j][str(values[i,locol+j])] = True

  grid = np.array([len(tempGrid[0]),\
                   len(tempGrid[1]),\
                   len(tempGrid[2]),],'i')
  
  dim = np.ones(3)

  for i,r in enumerate(grid):
    if r > 1:
      dim[i] = (max(map(float,tempGrid[i].keys()))-min(map(float,tempGrid[i].keys())))*r/(r-1.0)
  if grid[2]==1: # for 2D case set undefined dimension to given unitlength or alternatively give it the length of the smallest element
    if options.unitlength == 0.0: 
      dim[2] = min(dim/grid)
    else:
      dim[2] = options.unitlength
  print dim
  if options.undeformed:
    Favg = np.eye(3)
  else:
    Favg = damask.core.math.tensorAvg(
                      np.reshape(np.transpose(values[:,column['tensor'][options.defgrad]:
                                                             column['tensor'][options.defgrad]+9]),
                                                             (3,3,grid[0],grid[1],grid[2])))

  F = np.reshape(np.transpose(values[:,column['tensor'][options.defgrad]:
                                       column['tensor'][options.defgrad]+9]),
                                                           (3,3,grid[0],grid[1],grid[2]))
  centroids = damask.core.mesh.deformedCoordsFFT(dim,F,Favg,options.scaling)
  nodes = damask.core.mesh.nodesAroundCentres(dim,Favg,centroids)

  fields =  {\
             'tensor': {},\
             'vector': {},\
             'scalar': {},\
             'double': {},\
             'triple': {},\
             'quadruple': {},\
            }
  reshape = {\
             'tensor': [3,3],\
             'vector': [3],\
             'scalar': [],\
             'double': [2],\
             'triple': [3],\
             'quadruple': [4],\
            }
  length =  {\
             'tensor': 9,\
             'vector': 3,\
             'scalar': 1,\
             'double': 2,\
             'triple': 3,\
             'quadruple': 4,\
            }


# vtk lib out
  if False:
    points = vtk.vtkPoints()
    for z in range (grid[2]+1):
      for y in range (grid[1]+1):
        for x in range (grid[0]+1):
          points.InsertNextPoint(nodes[:,x,y,z])
   
    data=[]
    j=0
    for datatype in fields.keys():
      for what in eval('options.'+datatype):
        for label in matches[datatype][what]:
          col = column[datatype][label]
          if col != -1:
            data.append(vtk.vtkFloatArray())
            data[j].SetNumberOfComponents(length[datatype])
            for i in xrange(grid[2]*grid[1]*grid[0]):
              for k in xrange(length[datatype]):
                data[j].InsertNextValue(values2[i,col+k])
            data[j].SetName(label)
            j+=1
   
    if options.output_mesh:
      hexs = vtk.vtkCellArray()
      i = 0
      elems=[]
      for z in range (grid[2]):
        for y in range (grid[1]):
          for x in range (grid[0]):
   
            elems.append(vtk.vtkHexahedron())
            base = z*(grid[1]+1)*(grid[0]+1)+y*(grid[0]+1)+x
            elems[i].GetPointIds().SetId(0, base)
            elems[i].GetPointIds().SetId(1, base+1)
            elems[i].GetPointIds().SetId(2, base+grid[0]+2)
            elems[i].GetPointIds().SetId(3, base+grid[0]+1)
            elems[i].GetPointIds().SetId(4, base+(grid[1]+1)*(grid[0]+1))
            elems[i].GetPointIds().SetId(5, base+(grid[1]+1)*(grid[0]+1)+1)
            elems[i].GetPointIds().SetId(6, base+(grid[1]+1)*(grid[0]+1)+grid[0]+2)
            elems[i].GetPointIds().SetId(7, base+(grid[1]+1)*(grid[0]+1)+grid[0]+1)
            hexs.InsertNextCell(elems[i])
            i+=1
   
      uGrid = vtk.vtkUnstructuredGrid()
      uGrid.SetPoints(points)
      i = 0
      for z in range (grid[2]):
        for y in range (grid[1]):
          for x in range (grid[0]):
            uGrid.InsertNextCell(elems[i].GetCellType(), elems[i].GetPointIds())
            i+=1
   
      for i in xrange(len(data)):
        uGrid.GetCellData().AddArray(data[i])
   
      outWriter = vtk.vtkXMLUnstructuredGridWriter()
      outWriter.SetDataModeToBinary()
      outWriter.SetCompressorTypeToZLib()
      (head,tail) = os.path.split(filename)
      outWriter.SetFileName(os.path.join(head,'mesh_'+os.path.splitext(tail)[0]+'.vtu'))
      outWriter.SetInput(uGrid)
      outWriter.Write()

  
  for datatype in fields.keys():
    print '\n%s:'%datatype,
    fields[datatype]['_order_'] = []
    for what in eval('options.'+datatype):
      for label in matches[datatype][what]:
        col = column[datatype][label]
        if col != -1:
          print label,
          fields[datatype][label] = np.reshape(values[:,col:col+length[datatype]],[grid[0],grid[1],grid[2]]+reshape[datatype])
          fields[datatype]['_order_'] += [label]
  print '\n'

  out = {}
  if options.output_mesh:   out['mesh']   = vtk_writeASCII_mesh(nodes,fields,grid,sep[options.separator])
  if options.output_points: out['points'] = vtk_writeASCII_points(centroids,fields,grid,sep[options.separator])
  
  for what in out.keys():
    print what
    (head,tail) = os.path.split(filename)
    vtk = open(os.path.join(head,what+'_'+os.path.splitext(tail)[0]+'.vtk'), 'w')
    output(out[what],{'filepointer':vtk},'File')
    vtk.close()
  print
