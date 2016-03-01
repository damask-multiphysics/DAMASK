#!/usr/bin/env python

import os,sys,time,copy
import numpy as np
import damask
from optparse import OptionParser
from scipy import spatial
from collections import defaultdict

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
                  help = 'disorientation threshold per grain [%default] (degrees)')
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

parser.set_defaults(symmetry = 'cubic',
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
toRadians = np.pi/180.0 if options.degrees else 1.0                                               # rescale degrees to radians
cos_disorientation = np.cos(options.disorientation/2.*toRadians)

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header -------------------------------------------  

  table.head_read()

# ------------------------------------------ sanity checks -----------------------------------------

  errors  = []
  remarks = []
  
  if table.label_dimension(options.coords) != 3:       errors.append('coordinates {} are not a vector.'.format(options.coords))
  if not np.all(table.label_dimension(label) == dim):  errors.append('input {} has wrong dimension {}.'.format(label,dim))
  else:  column = table.label_index(label)

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append('grainID_{}@{}'.format(label,
                                             options.disorientation if options.degrees else np.degrees(options.disorientation)))    # report orientation source and disorientation in degrees
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
      bg.set_message('(%02i:%02i:%02i) processing point %i of %i (grain count %i)...'%(time_delta//3600,time_delta%3600//60,time_delta%60,p,len(grainID),len(orientations)))

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

    matched = False

# check against last matched needs to be really picky. best would be to exclude jumps across the poke (checking distance between last and me?)
# when walking through neighborhood first check whether grainID of that point has already been tested, if yes, skip!

    if matchedID != -1:                                                                             # has matched before?
      matched = (o.quaternion.conjugated() * orientations[matchedID].quaternion).w > cos_disorientation

    if not matched:
      alreadyChecked = {}
      bestDisorientation = damask.Quaternion([0,0,0,1])                                             # initialize to 180 deg rotation as worst case
      for i in kdtree.query_ball_point(kdtree.data[p],options.radius):                              # check all neighboring points
        gID = grainID[i]
        if gID != -1 and gID not in alreadyChecked:                                                 # an already indexed point belonging to a grain not yet tested?
          alreadyChecked[gID] = True                                                                # remember not to check again
          disorientation = o.disorientation(orientations[gID],SST = False)[0]                       # compare against that grain's orientation (and skip requirement of axis within SST)
          if disorientation.quaternion.w >  cos_disorientation and \
             disorientation.quaternion.w >= bestDisorientation.w:                                   # within disorientation threshold and better than current best?
            matched = True
            matchedID = gID                                                                         # remember that grain
            bestDisorientation = disorientation.quaternion

    if not matched:                                                                                 # no match -> new grain found
      memberCounts += [1]                                                                           # start new membership counter
      orientations += [o]                                                                           # initialize with current orientation
      matchedID = g
      g += 1                                                                                        # increment grain counter

    else:                                                                                           # did match existing grain
      memberCounts[matchedID] += 1

    grainID[p] = matchedID                                                                          # remember grain index assigned to point
    p += 1                                                                                          # increment point

  bg.set_message('identifying similar orientations among {} grains...'.format(len(orientations)))

  memberCounts = np.array(memberCounts)
  similarOrientations = [[] for i in xrange(len(orientations))]

  for i,orientation in enumerate(orientations[:-1]):                                                       # compare each identified orientation...
    for j in xrange(i+1,len(orientations)):                                                                # ...against all others that were defined afterwards
      if orientation.disorientation(orientations[j],SST = False)[0].quaternion.w > cos_disorientation:     # similar orientations in both grainIDs?
        similarOrientations[i].append(j)                                                                   # remember in upper triangle...
        similarOrientations[j].append(i)                                                                   # ...and lower triangle of matrix

    if similarOrientations[i] != []:
      bg.set_message('grainID {} is as: {}'.format(i,' '.join(map(str,similarOrientations[i]))))

  stillShifting = True
  while stillShifting:
    stillShifting = False
    tick = time.clock()

    for p,gID in enumerate(grainID):                                                                # walk through all points
      if p > 0 and p % 1000 == 0:

        time_delta = (time.clock()-tick) * (len(grainID) - p) / p
        bg.set_message('(%02i:%02i:%02i) shifting ID of point %i out of %i (grain count %i)...'%(time_delta//3600,time_delta%3600//60,time_delta%60,p,len(grainID),len(orientations)))
      if similarOrientations[gID] != []:                                                            # orientation of my grainID is similar to someone else?
        similarNeighbors = defaultdict(int)                                                         # dict holding frequency of neighboring grainIDs that share my orientation (freq info not used...)
        for i in kdtree.query_ball_point(kdtree.data[p],options.radius):                            # check all neighboring points
          if grainID[i] in similarOrientations[gID]:                                                # neighboring point shares my orientation?
            similarNeighbors[grainID[i]] += 1                                                       # remember its grainID
        if similarNeighbors != {}:                                                                  # found similar orientation(s) in neighborhood
          candidates = np.array([gID]+similarNeighbors.keys())                                      # possible replacement grainIDs for me
          grainID[p] = candidates[np.argsort(memberCounts[candidates])[-1]]                         # adopt ID that is most frequent in overall dataset
          memberCounts[gID]        -= 1                                                             # my former ID loses one fellow
          memberCounts[grainID[p]] += 1                                                             # my new ID gains one fellow
          bg.set_message('{}:{} --> {}'.format(p,gID,grainID[p]))                                   # report switch of grainID
          stillShifting = True

  table.data_rewind()

  outputAlive = True
  p = 0
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    table.data_append(1+grainID[p])                                                                 # add grain ID
    outputAlive = table.data_write()                                                                # output processed line
    p += 1

  bg.set_message('done after {} seconds'.format(time.clock()-start))

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables