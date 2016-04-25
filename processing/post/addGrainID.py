#!/usr/bin/env python

import os,sys,time,copy
import numpy as np
import damask
from optparse import OptionParser
from scipy import spatial

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add grain index based on similiarity of crystal lattice orientation.

""", version = scriptID)

parser.add_option('-r', '--radius',
                  dest = 'radius',
                  type = 'float', metavar = 'float',
                  help = 'search radius')
parser.add_option('-d', '--disorientation',
                  dest = 'disorientation',
                  type = 'float', metavar = 'float',
                  help = 'disorientation threshold in degrees [%default]')
parser.add_option('-s', '--symmetry',
                  dest = 'symmetry',
                  type = 'string', metavar = 'string',
                  help = 'crystal symmetry [%default]')
parser.add_option('-e', '--eulers',
                  dest = 'eulers',
                  type = 'string', metavar = 'string',
                  help = 'Euler angles')
parser.add_option(     '--degrees',
                  dest = 'degrees',
                  action = 'store_true',
                  help = 'Euler angles are given in degrees [%default]')
parser.add_option('-m', '--matrix',
                  dest = 'matrix',
                  type = 'string', metavar = 'string',
                  help = 'orientation matrix')
parser.add_option('-a',
                  dest = 'a',
                  type = 'string', metavar = 'string',
                  help = 'crystal frame a vector')
parser.add_option('-b',
                  dest = 'b',
                  type = 'string', metavar = 'string',
                  help = 'crystal frame b vector')
parser.add_option('-c',
                  dest = 'c',
                  type = 'string', metavar = 'string',
                  help = 'crystal frame c vector')
parser.add_option('-q', '--quaternion',
                  dest = 'quaternion',
                  type = 'string', metavar = 'string',
                  help = 'quaternion')
parser.add_option('-p', '--position',
                  dest = 'coords',
                  type = 'string', metavar = 'string',
                  help = 'spatial position of voxel [%default]')

parser.set_defaults(disorientation = 5,
                    symmetry = 'cubic',
                    coords   = 'pos',
                    degrees  = False,
                   )

(options, filenames) = parser.parse_args()

if options.radius is None:
  parser.error('no radius specified.')

input = [options.eulers     is not None,
         options.a          is not None and \
         options.b          is not None and \
         options.c          is not None,
         options.matrix     is not None,
         options.quaternion is not None,
        ]

if np.sum(input) != 1: parser.error('needs exactly one input format.')

(label,dim,inputtype) = [(options.eulers,3,'eulers'),
                         ([options.a,options.b,options.c],[3,3,3],'frame'),
                         (options.matrix,9,'matrix'),
                         (options.quaternion,4,'quaternion'),
                        ][np.where(input)[0][0]]                                                    # select input label that was requested
toRadians = np.pi/180.0 if options.degrees else 1.0                                                 # rescale degrees to radians
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
  
  if not 3 >= table.label_dimension(options.coords) >= 1:
    errors.append('coordinates "{}" need to have one, two, or three dimensions.'.format(options.coords))
  if not np.all(table.label_dimension(label) == dim):
    errors.append('input {} does not have dimension {}.'.format(label,dim))
  else:  column = table.label_index(label)

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append('grainID_{}@{:g}'.format('+'.join(label) 
                                               if isinstance(label, (list,tuple))
                                               else label,
                                               options.disorientation))                        # report orientation source and disorientation
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

# ------------------------------------------ build KD tree -----------------------------------------

# --- start background messaging

  bg = damask.util.backgroundMessage()
  bg.start()

  bg.set_message('reading positions...')

  table.data_readArray(options.coords)                                                              # read position vectors
  grainID = -np.ones(len(table.data),dtype=int)

  start = tick = time.clock()
  bg.set_message('building KD tree...')
  kdtree = spatial.KDTree(copy.deepcopy(table.data))

# ------------------------------------------ assign grain IDs --------------------------------------

  tick = time.clock()

  orientations = []                                                                                 # quaternions found for grain
  memberCounts = []                                                                                 # number of voxels in grain
  p = 0                                                                                             # point counter
  g = 0                                                                                             # grain counter
  matchedID = -1
  lastDistance = np.dot(kdtree.data[-1]-kdtree.data[0],kdtree.data[-1]-kdtree.data[0])              # (arbitrarily) use diagonal of cloud

  table.data_rewind()
  while table.data_read():                                                                          # read next data line of ASCII table

    if p > 0 and p % 1000 == 0:

      time_delta = (time.clock()-tick) * (len(grainID) - p) / p
      bg.set_message('(%02i:%02i:%02i) processing point %i of %i (grain count %i)...'\
            %(time_delta//3600,time_delta%3600//60,time_delta%60,p,len(grainID),np.count_nonzero(memberCounts)))

    if inputtype == 'eulers':
      o = damask.Orientation(Eulers   = np.array(map(float,table.data[column:column+3]))*toRadians,
                             symmetry = options.symmetry).reduced()
    elif inputtype == 'matrix':
      o = damask.Orientation(matrix   = np.array(map(float,table.data[column:column+9])).reshape(3,3).transpose(),
                             symmetry = options.symmetry).reduced()
    elif inputtype == 'frame':
      o = damask.Orientation(matrix = np.array(map(float,table.data[column[0]:column[0]+3] + \
                                                         table.data[column[1]:column[1]+3] + \
                                                         table.data[column[2]:column[2]+3])).reshape(3,3),
                             symmetry = options.symmetry).reduced()
    elif inputtype == 'quaternion':
      o = damask.Orientation(quaternion = np.array(map(float,table.data[column:column+4])),
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
        if disorientation.quaternion.w >  cos_disorientation:                                       # within threshold ...
          candidates.append(gID)                                                                    # remember as potential candidate
          if disorientation.quaternion.w >= bestDisorientation.w:                                   # ... and better than current best? 
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

  bg.set_message('done after {} seconds'.format(time.clock()-start))

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
