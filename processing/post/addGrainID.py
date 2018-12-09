#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-

import os,sys,time,copy
import numpy as np
import damask
from optparse import OptionParser
from scipy import spatial

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Add grain index based on similiarity of crystal lattice orientation.

""", version = scriptID)

parser.add_option('-r',
                  '--radius',
                  dest = 'radius',
                  type = 'float', metavar = 'float',
                  help = 'search radius')
parser.add_option('-d',
                  '--disorientation',
                  dest = 'disorientation',
                  type = 'float', metavar = 'float',
                  help = 'disorientation threshold in degrees [%default]')
parser.add_option('-s',
                  '--symmetry',
                  dest = 'symmetry',
                  type = 'string', metavar = 'string',
                  help = 'crystal symmetry [%default]')
parser.add_option('-q',
                  '--quaternion',
                  dest = 'quaternion',
                  type = 'string', metavar = 'string',
                  help = 'label of quaternion')
parser.add_option('-p',
                  '--pos', '--position',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'label of coordinates [%default]')

parser.set_defaults(disorientation = 5,
                    quaternion = 'orientation',
                    symmetry = 'cubic',
                    pos      = 'pos',
                   )

(options, filenames) = parser.parse_args()

if options.radius is None:
  parser.error('no radius specified.')

cos_disorientation = np.cos(np.radians(options.disorientation/2.))                                  # cos of half the disorientation angle

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header -------------------------------------------  

  table.head_read()

# ------------------------------------------ sanity checks -----------------------------------------

  errors  = []
  remarks = []
  
  if not 3 >= table.label_dimension(options.pos) >= 1:
    errors.append('coordinates "{}" need to have one, two, or three dimensions.'.format(options.pos))
  if not np.all(table.label_dimension(options.quaternion) == 4):
    errors.append('input "{}" does not have dimension 4.'.format(options.quaternion))
  else:  column = table.label_index(options.quaternion)

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append('grainID_{}@{:g}'.format('+'.join(options.quaternion),
                                               options.disorientation))                        # report orientation source and disorientation
  table.head_write()

# ------------------------------------------ build KD tree -----------------------------------------

  table.data_readArray(options.pos)                                                                 # read position vectors
  grainID = -np.ones(len(table.data),dtype=int)
  Npoints = table.data.shape[0]
  kdtree = spatial.KDTree(copy.deepcopy(table.data))

# ------------------------------------------ assign grain IDs --------------------------------------

  orientations = []                                                                                 # quaternions found for grain
  memberCounts = []                                                                                 # number of voxels in grain
  p = 0                                                                                             # point counter
  g = 0                                                                                             # grain counter
  matchedID = -1
  lastDistance = np.dot(kdtree.data[-1]-kdtree.data[0],kdtree.data[-1]-kdtree.data[0])              # (arbitrarily) use diagonal of cloud

  table.data_rewind()
  while table.data_read():                                                                          # read next data line of ASCII table

    if Npoints > 100 and p%(Npoints//100) == 0:                                                     # report in 1% steps if possible and avoid modulo by zero
      damask.util.print_progress(iteration=p,total=Npoints)

    o = damask.Orientation(quaternion = np.array(list(map(float,table.data[column:column+4]))),
                           symmetry   = options.symmetry).reduced()

    matched        = False
    alreadyChecked = {}
    candidates     = []
    bestDisorientation = damask.Quaternion([0,0,0,1])                                               # initialize to 180 deg rotation as worst case

    for i in kdtree.query_ball_point(kdtree.data[p],options.radius):                                # check all neighboring points
      gID = grainID[i]
      if gID != -1 and gID not in alreadyChecked:                                                   # indexed point belonging to a grain not yet tested?
        alreadyChecked[gID] = True                                                                  # remember not to check again
        disorientation = o.disorientation(orientations[gID],SST = False)[0]                         # compare against other orientation
        if disorientation.quaternion.q >  cos_disorientation:                                       # within threshold ...
          candidates.append(gID)                                                                    # remember as potential candidate
          if disorientation.quaternion.q >= bestDisorientation.q:                                   # ... and better than current best? 
            matched = True
            matchedID = gID                                                                         # remember that grain
            bestDisorientation = disorientation.quaternion

    if matched:                                                                                     # did match existing grain
      memberCounts[matchedID] += 1
      if len(candidates) > 1:                                                                       # ambiguity in grain identification?
        largestGrain = sorted(candidates,key=lambda x:memberCounts[x])[-1]                          # find largest among potential candidate grains
        matchedID    = largestGrain
        for c in [c for c in candidates if c != largestGrain]:                                      # loop over smaller candidates
          memberCounts[largestGrain] += memberCounts[c]                                             # reassign member count of smaller to largest
          memberCounts[c] = 0
        grainID = np.where(np.in1d(grainID,candidates), largestGrain, grainID)                      # relabel grid points of smaller candidates as largest one
      
    else:                                                                                           # no match -> new grain found
      orientations += [o]                                                                           # initialize with current orientation
      memberCounts += [1]                                                                           # start new membership counter
      matchedID = g
      g += 1                                                                                        # increment grain counter

    grainID[p] = matchedID                                                                          # remember grain index assigned to point
    p += 1                                                                                          # increment point

  grainIDs = np.where(np.array(memberCounts) > 0)[0]                                                # identify "live" grain identifiers
  packingMap = dict(zip(list(grainIDs),range(len(grainIDs))))                                       # map to condense into consecutive IDs
  
  table.data_rewind()

  outputAlive = True
  p = 0
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    table.data_append(1+packingMap[grainID[p]])                                                     # add (condensed) grain ID
    outputAlive = table.data_write()                                                                # output processed line
    p += 1

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
