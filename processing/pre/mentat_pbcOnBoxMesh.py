#!/usr/bin/env python

import sys,os,pwd,math,re,string, damask
from optparse import OptionParser

sys.path.append(damask.solver.Marc().libraryPath('../../'))

try:
  from py_mentat import *
except:
  print('error: no valid Mentat release found')
  sys.exit(-1)


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
      }[dest](cmd,locals)
  return



def servoLink():

  cmds = []
  base = ['x','y','z']
  box = {'min': {'x': float(sys.maxint),'y': float(sys.maxint),'z': float(sys.maxint)},
       'max': {'x':-float(sys.maxint),'y':-float(sys.maxint),'z':-float(sys.maxint)},
       'delta': {'x':0,'y':0,'z':0},
      }
  Nnodes = py_get_int("nnodes()")
  NodeCoords = [{'x':py_get_float("node_x(%i)"%(node)),
           'y':py_get_float("node_y(%i)"%(node)),
           'z':py_get_float("node_z(%i)"%(node)),} for node in range(1,1+Nnodes)]

  for node in range(Nnodes):    # find the bounding box
    for coord in base:          # check each direction in turn
      box['min'][coord] = min(box['min'][coord],NodeCoords[node][coord])
      box['max'][coord] = max(box['max'][coord],NodeCoords[node][coord])

  for coord in base:            # calc the dimension of the bounding box
    box['delta'][coord] = box['max'][coord] - box['min'][coord] 
  
  baseNode = {}
  linkNodes = []
  
  for node in range(Nnodes):    # loop over all nodes
    pos = {}
    key = {}
    maxFlag = {'x': False, 'y': False, 'z': False}
    Nmax = 0
    Nmin = 0
    for coord in base:                                # for each direction
      key[coord] = "%.8e"%NodeCoords[node][coord]     # translate position to string
      if (key[coord] == "%.8e"%box['min'][coord]):    # compare to min of bounding box (i.e. is on outer face?)
        Nmin += 1                                     # count outer (back) face membership
      elif (key[coord] == "%.8e"%box['max'][coord]):  # compare to max of bounding box (i.e. is on outer face?)
        Nmax += 1                                     # count outer (front) face memebership
        maxFlag[coord] = True                         # remember face membership (for linked nodes)

    if Nmin > 0 and Nmin > Nmax:                      # node is on more back than font faces
      # prepare for any non-existing entries in the data structure
      if key['x'] not in baseNode.keys():
        baseNode[key['x']] = {}
      if key['y'] not in baseNode[key['x']].keys():
        baseNode[key['x']][key['y']] = {}
      if key['z'] not in baseNode[key['x']][key['y']].keys():
        baseNode[key['x']][key['y']][key['z']] = 0
        
      baseNode[key['x']][key['y']][key['z']] = node+1   # remember the base node id

    elif Nmax > 0 and Nmax >= Nmin:                   # node is on at least as many front than back faces
      linkNodes.append({'id': node+1,'coord': NodeCoords[node], 'onFaces': Nmax,'faceMember': maxFlag})
  

  baseCorner = baseNode["%.8e"%box['min']['x']]["%.8e"%box['min']['y']]["%.8e"%box['min']['z']]     # detect ultimate base node
  
  for node in linkNodes:                        # loop over all linked nodes
    linkCoord = [node['coord']]                 # start list of control node coords with my coords
    for dir in base:                            # check for each direction
      if node['faceMember'][dir]:               # me on this front face
        linkCoord[0][dir] = box['min'][dir]     # project me onto rear face along dir
        linkCoord.append({'x':box['min']['x'],'y':box['min']['y'],'z':box['min']['z'],})      # append base corner
        linkCoord[-1][dir] = box['max'][dir]    # stretch it to corresponding control leg of "dir"

    nLinks = len(linkCoord)

    for dof in [1,2,3]:
      cmds.append([
        "*new_link *link_class servo",
        "*link_class servo *tied_node %i"%node['id'],
        "*link_class servo *tied_dof %i"%dof,
        "*servo_nterms %i"%(1+nLinks),
        ])
      for i in range(nLinks):
        cmds.append([
        "*link_class servo *servo_ret_node %i %i"%(i+1,baseNode["%.8e"%linkCoord[i]['x']]["%.8e"%linkCoord[i]['y']]["%.8e"%linkCoord[i]['z']]),
        "*link_class servo *servo_ret_dof %i %i"%(i+1,dof),
        "*link_class servo *servo_ret_coef %i 1"%(i+1),
        ])
      cmds.append([
      "*link_class servo *servo_ret_node %i %i"%(1+nLinks,baseCorner),
      "*link_class servo *servo_ret_dof %i %i"%(1+nLinks,dof),
      "*link_class servo *servo_ret_coef %i -%i"%(1+nLinks,nLinks-1),
      ])
  
  cmds.append([
    "*select_nodes",
    ["%i"%node['id'] for node in linkNodes],
    "#",
  ])
  
  return cmds

  

# ----------------------- MAIN -------------------------------
  
parser = OptionParser(usage='%prog [options]', description = """
Set up servo linking to achieve periodic boundary conditions for a regular hexahedral mesh presently opened in MSC.Mentat

""" + string.replace('$Id$','\n','\\n')
)

parser.add_option("-p", "--port", type="int",\
                                  dest="port",\
                                  help="Mentat connection port [%default]")
parser.add_option("-v", "--verbose", action="store_true",\
                                  dest="verbose",\
                                  help="write Mentat command stream also to stdout [%default]")
parser.set_defaults(port = 40007)
parser.set_defaults(verbose = False)

(options, args) = parser.parse_args()

outputLocals = {}
print 'waiting to connect...'
py_connect('',options.port)
output([\
        '*remove_all_servos',
        '*sweep_all',
        '*renumber_nodes',
        '*set_links off',
        ],outputLocals,'Mentat')     # script depends on consecutive numbering of nodes
cmds = servoLink()
print 'connected...'
output(cmds,outputLocals,'Mentat')
output([\
        '*set_links on',
        '*draw',
        ],outputLocals,'Mentat')     # script depends on consecutive numbering of nodes
py_disconnect()

if options.verbose:
  output(cmds,outputLocals,'Stdout')
