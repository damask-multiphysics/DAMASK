#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import sys,os,pwd,math,string
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

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
    sys.stdout.write(cmd+'\n')
  else:
    sys.stdout.write(cmd+'\n')
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
  box = {'min': np.zeros(3,dtype='d'),
         'max': np.zeros(3,dtype='d'),
       'delta': np.zeros(3,dtype='d'),
      }
  Nnodes = py_get_int("nnodes()")
  NodeCoords = np.zeros((Nnodes,3),dtype='d')
  for node in xrange(Nnodes):
    NodeCoords[node,0] = py_get_float("node_x(%i)"%(node+1))
    NodeCoords[node,1] = py_get_float("node_y(%i)"%(node+1))
    NodeCoords[node,2] = py_get_float("node_z(%i)"%(node+1))

  box['min'] = NodeCoords.min(axis=0)                   # find the bounding box
  box['max'] = NodeCoords.max(axis=0)
  box['delta'] = box['max']-box['min']
  for coord in xrange(3):                               # calc the dimension of the bounding box
    if box['delta'][coord] != 0.0:
      for extremum in ['min','max']:
        rounded = round(box[extremum][coord]*1e+15/box['delta'][coord]) * \
                                             1e-15*box['delta'][coord]       # rounding to 1e-15 of dimension
        box[extremum][coord] = {False: rounded,
                                 True: 0.0}[rounded == 0.0]                  # get rid of -0.0 (negative zeros)
  baseNode = {}
  linkNodes = []
  
  for node in xrange(Nnodes):                           # loop over all nodes
    pos = {}
    key = {}
    maxFlag = [False, False, False]
    Nmax = 0
    Nmin = 0
    for coord in xrange(3):                             # for each direction
      if box['delta'][coord] != 0.0:
        rounded = round(NodeCoords[node,coord]*1e+15/box['delta'][coord]) * \
                                               1e-15*box['delta'][coord]     # rounding to 1e-15 of dimension
        NodeCoords[node,coord] = {False: rounded,
                                   True: 0.0}[rounded == 0.0]                # get rid of -0.0 (negative zeros)
      key[base[coord]] = "%.8e"%NodeCoords[node,coord]                       # translate position to string
      if   (key[base[coord]] == "%.8e"%box['min'][coord]):                   # compare to min of bounding box (i.e. is on outer face?)
        Nmin += 1                                                            # count outer (back) face membership
      elif (key[base[coord]] == "%.8e"%box['max'][coord]):                   # compare to max of bounding box (i.e. is on outer face?)
        Nmax += 1                                                            # count outer (front) face membership
        maxFlag[coord] = True                                                # remember face membership (for linked nodes)

    if Nmin > 0 and Nmin > Nmax:                                             # node is on more back than front faces
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
  

  baseCorner = baseNode["%.8e"%box['min'][0]]["%.8e"%box['min'][1]]["%.8e"%box['min'][2]]     # detect ultimate base node
  
  
  for node in linkNodes:                          # loop over all linked nodes
    linkCoord = [node['coord']]                   # start list of control node coords with my coords
    for dir in xrange(3):                         # check for each direction
      if node['faceMember'][dir]:                 # me on this front face
        linkCoord[0][dir] = box['min'][dir]       # project me onto rear face along dir
        linkCoord.append(np.array(box['min'])) # append base corner
        
        linkCoord[-1][dir] = box['max'][dir]      # stretch it to corresponding control leg of "dir"

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
        "*link_class servo *servo_ret_node %i %i"%(i+1,baseNode["%.8e"%linkCoord[i][0]]["%.8e"%linkCoord[i][1]]["%.8e"%linkCoord[i][2]]),
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

#--------------------------------------------------------------------------------------------------
#                                MAIN
#-------------------------------------------------------------------------------------------------- 
parser = OptionParser(option_class=damask.extendableOption, usage = '%prog [options]', description = """
Set up servo linking to achieve periodic boundary conditions for a regular hexahedral mesh presently opened in MSC.Mentat

""", version = scriptID)

parser.add_option("-p", "--port", type="int", dest="port", metavar='int',
                                  help="Mentat connection port [%default]")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                                  help="write Mentat command stream also to stdout [%default]")
parser.set_defaults(port = 40007)
parser.set_defaults(verbose = False)

(options, args) = parser.parse_args()

outputLocals = {}
if options.verbose:
  file={'croak':sys.stderr}
else:
  file={'croak':sys.stdout}

file['croak'].write('\033[1m'+scriptName+'\033[0m\n\n')
file['croak'].write( 'waiting to connect...\n')
py_connect('',options.port)
file['croak'].write('connected...\n')
output([\
        '*draw_manual',              # prevent redrawing in "new" Mentat, should be much faster
        '*remove_all_servos',
        '*sweep_all',
        '*renumber_nodes',
        '*set_links off',
        ],outputLocals,'Mentat')     # script depends on consecutive numbering of nodes
cmds = servoLink()
output(cmds,outputLocals,'Mentat')
output([\
        '*set_links on',
        '*draw',
        ],outputLocals,'Mentat')     # script depends on consecutive numbering of nodes
py_disconnect()

if options.verbose:
  output(cmds,outputLocals,'Stdout')
