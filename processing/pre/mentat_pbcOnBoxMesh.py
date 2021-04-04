#!/usr/bin/env python3

import sys
import os
import re
import time
import tempfile
from optparse import OptionParser

import numpy as np

import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# Convert .mfd file into a usable format
# Broken into labeled sections (eg. nodes, links, etc)
# Each section has a list of labeled elements with formatted numerical data
def parseMFD(dat):
  formatted = []
  section = 0
  formatted.append({'label': 'header', 'uid': -1, 'els': []})
  # in between =beg= and =end= part of file
  in_block = False
  for line in dat:
    if in_block: # currently in a section
      # lines that start with a space are numerical data
      if line[0] == ' ':
        formatted[section]['els'].append([])

        # grab numbers
        nums = re.split(r'\s+', line.strip())

        for num in nums:
          # floating point has format ' -x.xxxxxxxxxxxxe+yy'
          # scientific notation is used for float
          if (len(num) >= 4) and (num[-4] == 'e'):
            formatted[section]['els'][-1].append(float(num))
          else: # integer
            formatted[section]['els'][-1].append(int(num))
      else: # not numerical data, so it is a label for an element or section end
        if line[0] == '=' and re.search(r'=end=$', line) is not None: # End of section, avoiding regex if possible
          in_block = False
        else:
          formatted[section]['els'].append([])
          formatted[section]['els'][-1] = line

    else: # Not in a section, we are looking for a =beg= now
      search = re.search(r'=beg=\s+(\d+)\s\((.*?)\)', line)
      if search is not None: # found start of a new section
        section += 1
        in_block = True
        formatted.append({'label': search.group(2), 'uid': int(search.group(1)), 'els': []})
      else: # No =beg= found, probably in the header
        # Either header or somthing we didn't plan for - just save the line so it isn't lost
        if formatted[section]['uid'] > 0:
          section += 1
          formatted.append({'label': '', 'uid': -2, 'els': []}) # make dummy section to store unrecognized data
        formatted[section]['els'].append(line)

  return formatted

def asMFD(mfd_data):
  result = ''
  for section in mfd_data:
    if section['uid'] > 0:
      result += '=beg={0:5d} ({1})\n'.format(section['uid'], section['label'])
    for el in section['els']:
      if type(el) == str:
        result += el
      elif type(el) == list:
        for num in el:
          if type(num) == int:
            result += '{:20d}'.format(num)
          elif type(num) == float:
            result += '{:20.12e}'.format(num)
          else:
            print(f'WARNING: encountered unknown type: {type(el)}')
        result += '\n'
      else:
        print(f'WARNING: encountered unknown type: {type(el)}')
    if section['uid'] > 0:
      result += '=end=\n'
  return result.strip()


def add_servoLinks(mfd_data,active=[True,True,True]):  # directions on which to add PBC
  base = ['x','y','z']
  box = {'min': np.zeros(3,dtype='d'),
         'max': np.zeros(3,dtype='d'),
       'delta': np.zeros(3,dtype='d'),
      }

  mfd_dict = {}
  for i in range(len(mfd_data)):
    mfd_dict[mfd_data[i]['label']] = i

  NodeCoords = np.array(mfd_data[mfd_dict['nodes']]['els'][1::4])[:,1:4]
  Nnodes = NodeCoords.shape[0]

  box['min'] = NodeCoords.min(axis=0)                                                               # find the bounding box
  box['max'] = NodeCoords.max(axis=0)
  box['delta'] = box['max']-box['min']
  for coord in range(3):                                                                            # calc the dimension of the bounding box
    if box['delta'][coord] != 0.0:
      for extremum in ['min','max']:
        rounded = round(box[extremum][coord]*1e+15/box['delta'][coord]) * \
                                             1e-15*box['delta'][coord]                              # rounding to 1e-15 of dimension
        box[extremum][coord] = 0.0 if rounded == 0.0 else rounded                                   # get rid of -0.0 (negative zeros)
  baseNode = {}
  linkNodes = []

#-------------------------------------------------------------------------------------------------
# loop over all nodes
  for node in range(Nnodes):
    key = {}
    maxFlag = [False, False, False]
    Nmax = 0
    Nmin = 0
    for coord in range(3):                                                                          # for each direction
      if box['delta'][coord] != 0.0:
        rounded = round(NodeCoords[node,coord]*1e+15/box['delta'][coord]) * \
                                               1e-15*box['delta'][coord]                            # rounding to 1e-15 of dimension
        NodeCoords[node,coord] = 0.0 if rounded == 0.0 else rounded                                 # get rid of -0.0 (negative zeros)
      key[base[coord]] = "%.8e"%NodeCoords[node,coord]                                              # translate position to string
      if   (key[base[coord]] == "%.8e"%box['min'][coord]):                                          # compare to min of bounding box (i.e. is on outer face?)
        Nmin += 1                                                                                   # count outer (back) face membership
      elif (key[base[coord]] == "%.8e"%box['max'][coord]):                                          # compare to max of bounding box (i.e. is on outer face?)
        Nmax += 1                                                                                   # count outer (front) face membership
        maxFlag[coord] = True                                                                       # remember face membership (for linked nodes)

    if Nmin > 0:                                                                                    # node is on a back face
        # prepare for any non-existing entries in the data structure
      if key['x'] not in baseNode.keys():
        baseNode[key['x']] = {}
      if key['y'] not in baseNode[key['x']].keys():
        baseNode[key['x']][key['y']] = {}
      if key['z'] not in baseNode[key['x']][key['y']].keys():
        baseNode[key['x']][key['y']][key['z']] = 0

      baseNode[key['x']][key['y']][key['z']] = node+1                                               # remember the base node id

    if Nmax > 0 and Nmax >= Nmin:                                                                   # node is on at least as many front than back faces
      if any([maxFlag[i] and active[i] for i in range(3)]):
        linkNodes.append({'id': node+1,'coord': NodeCoords[node], 'faceMember': [maxFlag[i] and active[i] for i in range(3)]})

  mfd_data[mfd_dict['entities']]['els'][0][0] += len(linkNodes) * 3

  baseCorner = baseNode["%.8e"%box['min'][0]]["%.8e"%box['min'][1]]["%.8e"%box['min'][2]]           # detect ultimate base node

  links = {'uid': 1705, 'label': 'links', 'els': [[7,0],[9,0]]}
  linkID = 0
  for node in linkNodes:                                                                            # loop over all linked nodes
    linkCoord = [node['coord']]                                                                     # start list of control node coords with my coords
    for dir in range(3):                                                                            # check for each direction
      if node['faceMember'][dir]:                                                                   # me on this front face
        linkCoord[0][dir] = box['min'][dir]                                                         # project me onto rear face along dir
        linkCoord.append(np.array(box['min']))                                                      # append base corner
        linkCoord[-1][dir] = box['max'][dir]                                                        # stretch it to corresponding control leg of "dir"

    nLinks = len(linkCoord)
    for dof in [1,2,3]:
      tied_node = node['id']
      nterms = 1 + nLinks

      linkID += 1
      # Link header
      links['els'].append('link{0}\n'.format(linkID))
      links['els'].append([linkID, 1])
      links['els'].append([0])
      links['els'].append([0])
      links['els'].append([0, 0, 0, tied_node])

      # these need to be put in groups of four
      link_payload = [dof, 0, nterms]

      # Individual node contributions (node, dof, coef.)
      for i in range(nterms):
        if i == nLinks:
          link_payload.append(baseCorner)
        else:
          link_payload.append(baseNode["%.8e"%linkCoord[i][0]]["%.8e"%linkCoord[i][1]]["%.8e"%linkCoord[i][2]])
      for i in range(nterms):
          link_payload.append(dof)
      for i in range(nterms):
        if i == nLinks:
          link_payload.append(1.0 - nLinks)
        else:
          link_payload.append(1.0)

      # Needs to be formatted 4 data points per row, character width of 20, so 80 total
      for j in range(0, len(link_payload), 4):
        links['els'].append(link_payload[j:j+4])
      if j+4 < len(link_payload):
        links['els'].append(link_payload[j+4:])

  i = 0
  while i < len(mfd_data) and mfd_data[i]['uid'] < 1705: i += 1

  if mfd_data[i]['uid'] == 1705: del mfd_data[i]
  mfd_data.insert(i, links)


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------

parser = OptionParser(usage='%prog options [file[s]]', description = """
Set up servo linking to achieve periodic boundary conditions for a regular hexahedral mesh.
Use *py_connection to operate on model presently opened in MSC.Mentat.
""", version = scriptID)

parser.add_option('-p', '--port',
                  type = int, metavar = 'int', default = None,
                  help = 'Mentat connection port')
parser.add_option('-x',
                  action = 'store_false', default = True,
                  help = 'no PBC along x direction')
parser.add_option('-y',
                  action = 'store_false', default = True,
                  help = 'no PBC along y direction')
parser.add_option('-z',
                  action = 'store_false', default = True,
                  help = 'no PBC along z direction')

(options, filenames) = parser.parse_args()

remote = options.port is not None

if remote and filenames != []:
  parser.error('file can not be specified when port is given.')
if filenames == []: filenames = [None]

if remote:
  sys.path.append(str(damask.solver.Marc().library_path))
  import py_mentat

  print(scriptName+': waiting to connect...')
  filenames = [os.path.join(tempfile._get_default_tempdir(), next(tempfile._get_candidate_names()) + '.mfd')]
  try:
    py_mentat.py_connect('',options.port)
    py_mentat.py_send('*set_save_formatted on')
    py_mentat.py_send('*save_as_model "{}" yes'.format(filenames[0]))
    py_mentat.py_get_int("nnodes()")
  except py_mentat.InputError as err:
    print(f'{err}. Try Tools/Python/"Run as Separate Process" & "Initiate".')
    sys.exit(-1)
  print( 'connected...')

for name in filenames:
  while remote and not os.path.exists(name): time.sleep(0.5)
  with  open( name,'r') if name is not None else sys.stdin as fileIn:
    print(scriptName+': '+name)
    mfd = parseMFD(fileIn)

  add_servoLinks(mfd,[options.x,options.y,options.z])
  with open( name,'w') if name is not None else sys.stdout as fileOut:
    fileOut.write(asMFD(mfd))

if remote:
  py_mentat.py_send('*open_model "{}"'.format(filenames[0]))
