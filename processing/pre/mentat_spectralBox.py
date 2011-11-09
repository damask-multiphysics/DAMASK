#!/usr/bin/env python

import os, sys, math, re, threading, time, string, msc_tools
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP


# ----------------------- FUNCTIONS ----------------------------

def outMentat(cmd,locals):
  if cmd[0:3] == '(!)':
    exec(cmd[3:])
  elif cmd[0:3] == '(?)':
    cmd = eval(cmd[3:])
    py_send(cmd)
  else:
    py_send(cmd)
  return

def outStdout(cmd,locals):
  if cmd[0:3] == '(!)':
    exec(cmd[3:])
  elif cmd[0:3] == '(?)':
    cmd = eval(cmd[3:])
    print cmd
  else:
    print cmd
  return


def output(cmds,locals,dest):
  for cmd in cmds:
    if isinstance(cmd,list):
      output(cmd,locals,dest)
    else:
      {\
      'Mentat': outMentat,\
      'Stdout': outStdout,\
      }[dest](str(cmd),locals)
  return


  
#--------------------
def init():
#--------------------
    return ["*new_model yes",
      "*reset",
      "*select_clear",
      "*set_element_class hex8",
      "*set_nodes off",
      "*elements_solid",
      "*show_view 4",
      "*reset_view",
      "*view_perspective",
      "*redraw",
      ]


#--------------------
def mesh(r,d):
#--------------------
    return [
  "*add_nodes",
  "%f %f %f"%(0.0,0.0,0.0),
  "%f %f %f"%(0.0,0.0,d[2]),
  "%f %f %f"%(0.0,d[1],d[2]),
  "%f %f %f"%(0.0,d[1],0.0),
  "%f %f %f"%(-d[0],0.0,0.0),
  "%f %f %f"%(-d[0],0.0,d[2]),
  "%f %f %f"%(-d[0],d[1],d[2]),
  "%f %f %f"%(-d[0],d[1],0.0),
  "*add_elements",
  range(1,9),
  "*sub_divisions",
  "%i %i %i"%(r[2],r[1],r[0]),
  "*subdivide_elements",
  "all_existing",
  "*set_sweep_tolerance",
  "%f"%(float(min(d))/max(r)/2.0),
  "*sweep_all",
  "*renumber_all",
  "*set_move_scale_factor x -1",
  "*move_elements",
  "all_existing",
  "*flip_elements",
  "all_existing",
  "*fill_view",
    ]


#--------------------
def material():
#--------------------
  cmds = [\
  "*new_mater standard",
  "*mater_option general:state:solid",
  "*mater_option structural:type:hypo_elast",
  "*mater_name",
  "hypela2",
  "*add_mater_elements",
  "all_existing",
  "*geometry_type mech_three_solid",
#  "*geometry_option red_integ_capacity:on",         # see below: reduced integration with one IP gave trouble being always OUTDATED...
  "*add_geometry_elements",
  "all_existing",
  ]
  
  return cmds
  

#--------------------
def geometry():
#--------------------
  cmds = [\
  "*geometry_type mech_three_solid",
#  "*geometry_option red_integ_capacity:on",
  "*add_geometry_elements",
  "all_existing",
  "*element_type 7",         # we are NOT using reduced integration (type 117) but opt for /elementhomogeneous/ in the respective phase description (material.config)
  "all_existing",
  ]
  
  return cmds
  

#--------------------
def initial_conditions(homogenization,grains):
#--------------------
  elements = []
  element = 0
  for id in grains:
    element += 1 
    if len(elements) < id:
      for i in range(id-len(elements)):
        elements.append([])
    elements[id-1].append(element)

  cmds = [\
    "*new_icond",
    "*icond_name _temperature",
    "*icond_type state_variable",
    "*icond_param_value state_var_id 1",
    "*icond_dof_value var 300",
    "*add_icond_elements",
    "all_existing",
    "*new_icond",
    "*icond_name _homogenization",
    "*icond_type state_variable",
    "*icond_param_value state_var_id 2",
    "*icond_dof_value var %i"%homogenization,
    "*add_icond_elements",
    "all_existing",
    ]

  for grain,elementList in enumerate(elements):
    cmds.append([\
            "*new_icond",
            "*icond_name microstructure_%i"%(grain+1),
            "*icond_type state_variable",
            "*icond_param_value state_var_id 3",
            "*icond_dof_value var %i"%(grain+1),
            "*add_icond_elements",
            elementList,
            "#",
                ])
  return cmds


#--------------------
def parse_geomFile(content,homog):
#--------------------
  (skip,key) = content[0].split()[:2]
  if key[:4].lower() == 'head':
    skip = int(skip)+1
  else:
    skip = 0
  
  res = [0,0,0]
  dim = [0.0,0.0,0.0]
  homog = 0
  
  for line in content[:skip]:
    data = line.split()
    if data[0].lower() == 'resolution':
      res = map(int,data[2:8:2])
    if data[0].lower() == 'dimension':
      dim = map(float,data[2:8:2])
    if data[0].lower() == 'homogenization':
      homog = int(data[1])

  grains = []
  for line in content[skip:]:
    grains.append(int(line.split()[0]))

  return (res,dim,homog,grains)

#--------------------
def parse_spectralFile(content,homog):
#--------------------

  coords = [{},{},{}]
  maxBox = [-1.0e20,-1.0e20,-1.0e20]
  minBox = [ 1.0e20, 1.0e20, 1.0e20]
  dim = [0.0,0.0,0.0]
  res = [0,0,0]
  grains = []
  
  for line in content:
    data = line.split()[3:7]
    grains.append(int(data[3]))
    for i in range(3):
      maxBox[i] = max(maxBox[i],float(data[i]))
      minBox[i] = min(minBox[i],float(data[i]))
      coords[i][data[i]] = True
      
  for i in range(3):
    res[i] = len(coords[i])
    dim[i] = (maxBox[i]-minBox[i])*res[i]/(res[i]-1.0)

  return (res,dim,homog,grains)

# ----------------------- MAIN -------------------------------

parser = OptionParser(usage='%prog [options] spectral.datafile', description = """
Generate FE hexahedral mesh from spectral description file.

Acceptable formats are
geom: header plus list of grain numbers or
spectral: phi1,Phi,phi2,x,y,z,id,phase.

""" + string.replace('$Id$','\n','\\n')
)
parser.add_option("-p", "--port", type="int",\
                                  dest="port",\
                                  help="Mentat connection port")
parser.add_option("-g", "--geom", action="store_const", const="geom",\
                                  dest="filetype",\
                                  help="file has 'geom' format")
parser.add_option("-s", "--spectral", action="store_const", const="spectral",\
                                      dest="filetype",\
                                      help="file has 'spectral' format")

parser.add_option("--homogenization", type="int",\
                                      dest="homogenization",\
                                      help="homogenization index from material.config (only required for spectral file type)")


parser.set_defaults(filetype = 'geom')
parser.set_defaults(homogenization = 1)

(options, args) = parser.parse_args()


sys.path.append(msc_tools.MSC_TOOLS().library_paths(sys.argv[0],'../../'))

try:
  from py_mentat import *
except:
  print('no valid Mentat release found')
  if options.port != None: sys.exit(-1)

if not os.path.isfile(args[0]):
  parser.error("cannot open %s"%args[0])

file = open(args[0])
content = file.readlines()
file.close()

if options.filetype not in ['spectral','geom']:
  options.filetype = os.path.splitext(args[0])[1][1:]

print '\nparsing %s...'%options.filetype,
sys.stdout.flush()

(res,dim,homog,grains) = {\
  'geom':     parse_geomFile,
  'spectral': parse_spectralFile,
  }[options.filetype](content,options.homogenization)

print '%i grains in %s with resolution %s and homogenization %i\n'%(len(list(set(grains))),str(dim),str(res),homog)


cmds = [\
  init(),
  mesh(res,dim),
  material(),
  geometry(),
  initial_conditions(homog,grains),
  '*identify_sets',
  '*redraw',
]

outputLocals = {}
if (options.port != None):
  py_connect('',options.port)
  output(cmds,outputLocals,'Mentat')
  py_disconnect()
else:
  output(cmds,outputLocals,'Stdout')

