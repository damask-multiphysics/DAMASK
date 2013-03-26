#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

# This script is used for the post processing of the results achieved by the spectral method.
# As it reads in the data coming from "materialpoint_results", it can be adopted to the data
# computed using the FEM solvers. Its capable to handle elements with one IP in a regular order

import os,sys,threading,re,numpy,time,string,damask
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP

# -----------------------------
class extendedOption(Option):
# -----------------------------
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

            
# -----------------------------
class backgroundMessage(threading.Thread):
# -----------------------------
    
    def __init__(self):
        threading.Thread.__init__(self)
        self.message = ''
        self.new_message = ''
        self.counter = 0
        self.symbols = ['- ', '\ ', '| ', '/ ']
        self.waittime = 0.5
    
    def __quit__(self):
        length = len(self.message) + len(self.symbols[self.counter])
        sys.stderr.write(chr(8)*length + ' '*length + chr(8)*length)
        sys.stderr.write('')
    
    def run(self):
        while not threading.enumerate()[0]._Thread__stopped:
            time.sleep(self.waittime)
            self.update_message()
        self.__quit__()

    def set_message(self, new_message):
        self.new_message = new_message
        self.print_message()
    
    def print_message(self):
        length = len(self.message) + len(self.symbols[self.counter])
        sys.stderr.write(chr(8)*length + ' '*length + chr(8)*length)                                # delete former message
        sys.stderr.write(self.symbols[self.counter] + self.new_message)                             # print new message
        self.message = self.new_message
        
    def update_message(self):
        self.counter = (self.counter + 1)%len(self.symbols)
        self.print_message()



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
          string.replace('powered by $Id$','\n','\\n'),
          'ASCII',
          'DATASET UNSTRUCTURED_GRID',
          'POINTS %i float'%N1,
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
    for item in data[type]:
      cmds += [\
               '%s %s float %i'%(info[type]['name'].upper()+plural,item,info[type]['len']),
               {True:'LOOKUP_TABLE default',False:''}[info[type]['name'][:3]=='sca'],
               [[[sep.join(map(unravel,data[type][item][:,j,k]))] for j in range(res[1])] for k in range(res[2])],
              ]

  return cmds
  
# ++++++++++++++++++++++++++++++++++++++++++++++++++++
def gmsh_writeASCII_mesh(mesh,data,res,sep):
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
          '$MeshFormat',
          '2.1 0 8',
          '$EndMeshFormat',
          '$Nodes',
          '%i float'%N1,
          [[['\t'.join(map(str,l,mesh[:,i,j,k])) for l in range(1,N1+1) for i in range(res[0]+1)] for j in range(res[1]+1)] for k in range(res[2]+1)],
          '$EndNodes',
          '$Elements',
          '%i'%N,
          ]

# cells
  n_elem = 0
  for z in range (res[2]):
    for y in range (res[1]):
      for x in range (res[0]):
        base = z*(res[1]+1)*(res[0]+1)+y*(res[0]+1)+x
        n_elem +=1
        cmds.append('\t'.join(map(str,[ \
                                            n_elem,
                                            '5',
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
           'ElementData',
           '1',
           '%s'%item,     # name of the view
           '0.0',         # thats the time value
           '3', 
           '0',           # time step
           '1',
           '%i'%N
          ]
  
  for type in data:
    plural = {True:'',False:'S'}[type.lower().endswith('s')]
    for item in data[type]:
      cmds += [\
               '%s %s float %i'%(info[type]['name'].upper()+plural,item,info[type]['len']),
               'LOOKUP_TABLE default',
               [[[sep.join(map(str,data[type][item][:,j,k]))] for j in range(res[1])] for k in range(res[2])],
              ]

  return cmds
 
# +++++++++++++++++++++++++++++++++++++++++++++++++++
def vtk_writeASCII_points(coordinates,data,res,sep):
# +++++++++++++++++++++++++++++++++++++++++++++++++++
  """ function writes data array defined on a point field """
  N  = res[0]*res[1]*res[2]
  
  cmds = [\
          '# vtk DataFile Version 3.1',
          'powered by $Id$',
          'ASCII',
          'DATASET UNSTRUCTURED_GRID',
          'POINTS %i float'%N,
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
               '%s %s float'%(type.upper()+plural,item),
               {True:'LOOKUP_TABLE default',False:''}[type.lower()[:3]=='sca'],
               [[[sep.join(map(unravel,data[type][item][:,j,k]))] for j in range(res[1])] for k in range(res[2])],
              ]

  return cmds
  
# +++++++++++++++++++++++++++++++++++++++++++++++++++++
def vtk_writeASCII_box(diag,defgrad):
# +++++++++++++++++++++++++++++++++++++++++++++++++++++
  """ corner box for the average defgrad """
  points = numpy.array([\
                         [0.0,0.0,0.0,],\
                         [diag[0],0.0,0.0,],\
                         [diag[0],diag[1],0.0,],\
                         [0.0,diag[1],0.0,],\
                         [0.0,0.0,diag[2],],\
                         [diag[0],0.0,diag[2],],\
                         [diag[0],diag[1],diag[2],],\
                         [0.0,diag[1],diag[2],],\
                       ])

  cmds = [\
    '# vtk DataFile Version 3.1',
    'powered by $Id$',
    'ASCII',
    'DATASET UNSTRUCTURED_GRID',
    'POINTS 8 float',
    ['\t'.join(map(str,numpy.dot(defgrad_av,points[p]))) for p in range(8)],
    'CELLS 8 16',
    ['1\t%i'%i for i in range(8)],
    'CELL_TYPES 8',
    ['1']*8,
  ]
  
  return cmds



# ----------------------- MAIN -------------------------------

parser = OptionParser(option_class=extendedOption, usage='%prog [options] datafile[s]', description = """
Produce VTK file from data field. Coordinates are taken from (consecutive) x, y, and z columns.

""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-s', '--scalar', action='extend', dest='scalar', type='string', \
                  help='list of scalars to visualize')
parser.add_option(      '--double', action='extend', dest='double', type='string', \
                  help='list of scalars to visualize')
parser.add_option(      '--triple', action='extend', dest='triple', type='string', \
                  help='list of scalars to visualize')
parser.add_option(      '--quadruple', action='extend', dest='quadruple', type='string', \
                  help='list of scalars to visualize')
parser.add_option('-v', '--vector', action='extend', dest='vector', type='string', \
                  help='list of vectors to visualize')
parser.add_option('-d', '--deformation', dest='defgrad', type='string', \
                  help='heading of deformation gradient columns [%default]')
parser.add_option('--reference', dest='undeformed', action='store_true',\
                  help='map results to reference (undeformed) configuration')
parser.add_option('-c','--cell', dest='cell', action='store_true',\
                  help='data is cell-centered')
parser.add_option('-p','--vertex', dest='cell', action='store_false',\
                  help='data is vertex-centered')
parser.add_option('--mesh', dest='output_mesh', action='store_true', \
                  help='produce VTK mesh file [%default]')
parser.add_option('--nomesh', dest='output_mesh', action='store_false', \
                  help='omit VTK mesh file')
parser.add_option('--points', dest='output_points', action='store_true', \
                  help='produce VTK points file [%default]')
parser.add_option('--nopoints', dest='output_points', action='store_false', \
                  help='omit VTK points file')
parser.add_option('--box', dest='output_box', action='store_true', \
                  help='produce VTK box file [%default]')
parser.add_option('--nobox', dest='output_box', action='store_false', \
                  help='omit VTK box file')
parser.add_option('--separator', dest='separator', type='string', \
                  help='data separator [t(ab), n(ewline), s(pace)]')
parser.add_option('--scaling', dest='scaling', action='extend', type='string', \
                  help='scaling of fluctuation')
parser.add_option('-u', '--unitlength', dest='unitlength', type='float', \
                  help='set unit length for 2D model [%default]')
parser.add_option('--filenodalcoords', dest='filenodalcoords', type='string', \
                  help='ASCII table containing nodal coords')
parser.add_option('--labelnodalcoords', dest='nodalcoords', type='string', nargs=3, \
                  help='labels of nodal coords in ASCII table')
parser.add_option('-l', '--linear', dest='linearreconstruction', action='store_true',\
                  help='use linear reconstruction of geometry [%default]')
                  
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
parser.set_defaults(output_box = False)
parser.set_defaults(scaling = [])
parser.set_defaults(undeformed = False)
parser.set_defaults(unitlength = 0.0)
parser.set_defaults(cell = True)
parser.set_defaults(filenodalcoords = '')
parser.set_defaults(labelnodalcoords = ['coord.x','coord.y','coord.z'])
parser.set_defaults(linearreconstruction = False)

sep = {'n': '\n', 't': '\t', 's': ' '}

(options, args) = parser.parse_args()

options.scaling += numpy.ones(max(0,3-len(options.scaling)))
options.scaling = map(float, options.scaling)
if numpy.any(options.scaling != 1.0) and options.linearreconstruction:  print 'cannot scale for linear reconstruction'
if numpy.any(options.scaling != 1.0) and options.filenodalcoords != '': print 'cannot scale when reading coordinate from file'
options.separator = options.separator.lower()
for filename in args:
  if not os.path.exists(filename):
    continue
  file = open(filename)
  content = file.readlines()
  file.close()
  m = re.search('(\d+)\shead',content[0],re.I)
  if m == None:
    continue
  print filename,'\n'
  sys.stdout.flush()
  
  headrow = int(m.group(1))
  headings = content[headrow].split()
  column = {}
  maxcol = 0
  
  for col,head in enumerate(headings):
    if head == {True:'ip.x',False:'node.x'}[options.cell]:
      locol = col
      maxcol = max(maxcol,col+3)
      break

  if locol < 0:
    print 'missing coordinates..!'
    continue
    
  column['tensor'] = {}
  for label in [options.defgrad] + options.tensor:
    column['tensor'][label] = -1
    for col,head in enumerate(headings):
      if head == label or head == '1_'+label:
        column['tensor'][label] = col
        maxcol = max(maxcol,col+9)
        break
      
  if not options.undeformed and column['tensor'][options.defgrad] < 0:
    print 'missing deformation gradient "%s"..!'%options.defgrad
    continue

  column['vector'] = {}
  for label in options.vector:
    column['vector'][label] = -1
    for col,head in enumerate(headings):
      if head == label or head == '1_'+label:
        column['vector'][label] = col
        maxcol = max(maxcol,col+3)
        break

  for length,what in enumerate(['scalar','double','triple','quadruple']):
    column[what] = {}
    for label in eval('options.'+what):
      column[what][label] = -1
      for col,head in enumerate(headings):
        if head == label or head == '1_'+label:
          column[what][label] = col
          maxcol = max(maxcol,col+1+length)
          break


  values = numpy.array(sorted([map(transliterateToFloat,line.split()[:maxcol]) for line in content[headrow+1:]],
                              key=lambda x:(x[locol+0],x[locol+1],x[locol+2])),'d')             # sort with z as fastest and x as slowest index

  N = len(values)

  grid = [{},{},{}]
  for j in xrange(3):
    for i in xrange(N):
      grid[j][str(values[i,locol+j])] = True

  res = numpy.array([len(grid[0]),\
                     len(grid[1]),\
                     len(grid[2]),],'i')
  if (res[0]%2 != 0 or res[1]%2 != 0 or (res[2] != 1 and res[2]%2 !=0)):
    print 'using linear reconstruction for uneven resolution'
    options.linearreconstruction = True
  dim = numpy.ones(3)

  for i,r in enumerate(res):
    if r > 1:
      dim[i] = (max(map(float,grid[i].keys()))-min(map(float,grid[i].keys())))*r/(r-1.0)
  if res[2]==1: # for 2D case set undefined dimension to given unitlength or alternatively give it the length of the smallest element
    if options.unitlength == 0.0: 
      dim[2] = min(dim/res)
    else:
      dim[2] = options.unitlength
 
  if options.undeformed:
    Favg = numpy.eye(3)
  else:
    Favg = damask.core.math.tensorAvg(
                      numpy.reshape(numpy.transpose(values[:,column['tensor'][options.defgrad]:
                                                             column['tensor'][options.defgrad]+9]),
                                                             (3,3,res[0],res[1],res[2])))
  if not options.filenodalcoords:
    F = numpy.reshape(numpy.transpose(values[:,column['tensor'][options.defgrad]:
                                               column['tensor'][options.defgrad]+9]),
                                                             (3,3,res[0],res[1],res[2]))
    if options.linearreconstruction:
      centroids = damask.core.mesh.deformedCoordsLinear(dim,F,Favg)
    else:
      centroids = damask.core.mesh.deformedCoordsFFT(dim,F,Favg,options.scaling)
    mesh = damask.core.mesh.nodesAroundCentres(dim,Favg,centroids)

  else:
    mesh = numpy.zeros(((res[0]+1)*(res[1]+1)*(res[2]+1),3),'d')

    filenodalcoords = open(options.filenodalcoords)
    tablenodalcoords = damask.ASCIItable(filenodalcoords)
    tablenodalcoords.head_read()
    columns = [tablenodalcoords.labels.index(options.labelnodalcoords[0]),        
               tablenodalcoords.labels.index(options.labelnodalcoords[1]),
               tablenodalcoords.labels.index(options.labelnodalcoords[2])]
    i = 0
    while tablenodalcoords.data_read():                                                    # read next data line of ASCII table


      mesh[i,:]=float(tablenodalcoords.data[column[:]])


      i += 1
    mesh=mesh.reshape(res[0]+1,res[1]+1,res[2]+1,3)

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

  for datatype in fields.keys():
    print '\n%s:'%datatype,
    for what in eval('options.'+datatype):
      col = column[datatype][what]
      if col != -1:
        print what,
        fields[datatype][what] = numpy.reshape(values[:,col:col+length[datatype]],[res[0],res[1],res[2]]+reshape[datatype])
  print '\n'

  out = {}
  if options.output_mesh:   out['mesh']   = vtk_writeASCII_mesh(mesh,fields,res,sep[options.separator])
  if options.output_points: out['points'] = vtk_writeASCII_points(centroids,fields,res,sep[options.separator])
  if options.output_box:    out['box']    = vtk_writeASCII_box(dim,Favg)
  
  for what in out.keys():
    print what
    (head,tail) = os.path.split(filename)
    vtk = open(os.path.join(head,what+'_'+os.path.splitext(tail)[0]+'.vtk'), 'w')
    output(out[what],{'filepointer':vtk},'File')
    vtk.close()
  print
  
