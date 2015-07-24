#!/usr/bin/env python

import os,sys,string,itertools,re,time,copy,operator,threading
import numpy as np
import damask
from scipy import spatial
from collections import defaultdict
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP

scriptID   = string.replace('$Id: addGrainID.py 2549 2013-07-10 09:13:21Z MPIE\p.eisenlohr $','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

#--------------------------------------------------------------------------------------------------
class extendedOption(Option):
#--------------------------------------------------------------------------------------------------
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
    self.symbols = ['- ', '\ ', '| ', '/ ',]
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
    sys.stderr.write(chr(8)*length + ' '*length + chr(8)*length)                  # delete former message
    sys.stderr.write(self.symbols[self.counter] + self.new_message)               # print new message
    self.message = self.new_message

  def update_message(self):
    self.counter = (self.counter + 1)%len(self.symbols)
    self.print_message()


parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
Add grain index based on similiarity of crystal lattice orientation.
""" + string.replace(scriptID,'\n','\\n')
)

parser.add_option('-r', '--radius', dest='radius', type='float',
                  help = 'search radius')
parser.add_option('-d', '--disorientation', dest='disorientation', type='float', metavar='ANGLE',
                  help = 'disorientation threshold per grain [%default] (degrees)')
parser.add_option('-s', '--symmetry', dest='symmetry', type='string',
                  help = 'crystal symmetry [%default]')
parser.add_option('-e', '--eulers',   dest='eulers', type='string', metavar='LABEL',
                  help = 'Euler angles')
parser.add_option(     '--degrees',   dest='degrees', action='store_true',
                  help = 'Euler angles are given in degrees [%default]')
parser.add_option('-m', '--matrix',   dest='matrix', type='string', metavar='LABEL',
                  help = 'orientation matrix')
parser.add_option('-a',               dest='a', type='string', metavar='LABEL',
                  help = 'crystal frame a vector')
parser.add_option('-b',               dest='b', type='string', metavar='LABEL',
                  help = 'crystal frame b vector')
parser.add_option('-c',               dest='c', type='string', metavar='LABEL',
                  help = 'crystal frame c vector')
parser.add_option('-q', '--quaternion', dest='quaternion', type='string', metavar='LABEL',
                  help = 'quaternion')
parser.add_option('-p', '--position', dest='position', type='string', metavar='LABEL',
                  help = 'spatial position of voxel [%default]')

parser.set_defaults(symmetry = 'cubic')
parser.set_defaults(position = 'pos')
parser.set_defaults(degrees = False)

(options, filenames) = parser.parse_args()

if options.radius == None:
  parser.error('no radius specified.')

datainfo = {                                                               # list of requested labels per datatype
             'tensor':     {'len':9,
                            'label':[]},
             'vector':     {'len':3,
                            'label':[]},
             'quaternion': {'len':4,
                            'label':[]},
           }

if options.eulers     != None:  datainfo['vector']['label'] += [options.eulers];                input = 'eulers'
if options.a          != None and \
   options.b          != None and \
   options.c          != None:  datainfo['vector']['label'] += [options.a,options.b,options.c]; input = 'frame'
if options.matrix     != None:  datainfo['tensor']['label'] += [options.matrix];                input = 'matrix'
if options.quaternion != None:  datainfo['quaternion']['label'] += [options.quaternion];        input = 'quaternion'

datainfo['vector']['label'] += [options.position]

toRadians = np.pi/180.0 if options.degrees else 1.0                                                                               # rescale degrees to radians
cos_disorientation = np.cos(options.disorientation/2.0*toRadians)

# ------------------------------------------ setup file handles ---------------------------------------

files = []
if filenames == []:
  files.append({'name':'STDIN',
                'input':sys.stdin,
                'output':sys.stdout,
                'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name,
                    'input':open(name),
                    'output':open(name+'_tmp','w'),
                    'croak':sys.stderr})

#--- loop over input files ------------------------------------------------------------------------

for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  table = damask.ASCIItable(file['input'],file['output'],buffered = False)  # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info

# --------------- figure out columns to process

  column = {}
  missingColumns = False

  for datatype,info in datainfo.items():
    for label in info['label']:
      key = list(set([label, '1_'+label]) & set(table.labels))                                      # check for intersection with table labels
      if key == []:
        file['croak'].write('column %s not found...\n'%label)
        missingColumns = True                                                                       # break if label not found
      else:
        column[label] = table.labels.index(key[0])                                                  # remember columns of requested data

  if missingColumns:
    continue

  table.labels_append('grainID_%g'%options.disorientation)

# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(string.replace(scriptID,'\n','\\n') + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

# ------------------------------------------ process data ---------------------------------------

# ------------------------------------------ build KD tree ---------------------------------------


# --- start background messaging

  bg = backgroundMessage()
  bg.start()

  bg.set_message('reading positions...')

  backup_readSize = table.__IO__['validReadSize']          # bad hack to circumvent overwriting by readArray...
  backup_labels = table.labels                             # bad hack...
  table.data_rewind()
  table.data_readArray(range(column[options.position],
                             column[options.position]+datainfo['vector']['len']))   # read position vectors
#  file['croak'].write('%i\n'%(len(table.data)))
  grainID = -np.ones(len(table.data),dtype=int)

  start = tick = time.clock()
  bg.set_message('building KD tree...')
  kdtree = spatial.KDTree(copy.deepcopy(table.data))
#  neighborhood = kdtree.query_ball_tree(kdtree,options.radius)
#  file['croak'].write('%.2f seconds\n'%(time.clock()-tick))
#  file['croak'].write('%i points\n'%(len(neighborhood)))


# ------------------------------------------ assign grain IDs ---------------------------------------

  orientations = []                                                         # quaternions found for grain
  memberCounts = []                                                         # number of voxels in grain

  table.data_rewind()
  table.__IO__['validReadSize'] = backup_readSize           # bad hack to circumvent overwriting by readArray...
  table.labels = backup_labels                              # bad hack...
  p = 0                                                                     # point counter
  g = 0                                                                     # grain counter
  matchedID = -1
  lastDistance = np.dot(kdtree.data[-1]-kdtree.data[0],kdtree.data[-1]-kdtree.data[0]) # (arbitrarily) use diagonal of cloud

  tick = time.clock()
  while table.data_read():                                                  # read next data line of ASCII table

    if p > 0 and p % 1000 == 0:

      time_delta = (time.clock()-tick) * (len(grainID) - p) / p
      bg.set_message('(%02i:%02i:%02i) processing point %i of %i (grain count %i)...'%(time_delta//3600,time_delta%3600//60,time_delta%60,p,len(grainID),len(orientations)))

    if input == 'eulers':
      o = damask.Orientation(Eulers=toRadians*\
                 np.array(map(float,table.data[column[options.eulers]:\
                                               column[options.eulers]+datainfo['vector']['len']])),
                             symmetry=options.symmetry).reduced()
    elif input == 'matrix':
      o = damask.Orientation(matrix=\
                 np.array([map(float,table.data[column[options.matrix]:\
                                                column[options.matrix]+datainfo['tensor']['len']])]).reshape(np.sqrt(datainfo['tensor']['len']),
                                                                                                             np.sqrt(datainfo['tensor']['len'])).transpose(),
                             symmetry=options.symmetry).reduced()
    elif input == 'frame':
      o = damask.Orientation(matrix=\
                 np.array([map(float,table.data[column[options.a]:\
                                                column[options.a]+datainfo['vector']['len']] + \
                                     table.data[column[options.b]:\
                                                column[options.b]+datainfo['vector']['len']] + \
                                     table.data[column[options.c]:\
                                                column[options.c]+datainfo['vector']['len']]
                                                    )]).reshape(3,3),
                             symmetry=options.symmetry).reduced()
    elif input == 'quaternion':
      o = damask.Orientation(quaternion=\
                 np.array(map(float,table.data[column[options.quaternion]:\
                                               column[options.quaternion]+datainfo['quaternion']['len']])),
                             symmetry=options.symmetry).reduced()

    matched = False

# check against last matched needs to be really picky. best would be to exclude jumps across the poke (checking distance between last and me?)
# when walking through neighborhood first check whether grainID of that point has already been tested, if yes, skip!

    if matchedID != -1:                                                       # has matched before?
      matched = (o.quaternion.conjugated() * orientations[matchedID].quaternion).w > cos_disorientation
#     if matchedID > 0:                                                       # has matched before?
#       thisDistance = np.dot(kdtree.data[p]-kdtree.data[p-1],kdtree.data[p]-kdtree.data[p-1],)
#       if thisDistance < 4.*lastDistance:                                    # about as close as last point pair?
#         disorientation = o.disorientation(orientations[matchedID-1]).quaternion.w  # check whether former grainID matches now again
#         matched = disorientation > cos_disorientation
#       lastDistance = thisDistance
#

    if not matched:
      alreadyChecked = {}
      bestDisorientation = damask.Orientation(quaternion=np.array([0,0,0,1]),symmetry = options.symmetry)  # initialize to 180 deg rotation as worst case
      for i in kdtree.query_ball_point(kdtree.data[p],options.radius):      # check all neighboring points
        gID = grainID[i]
        if gID != -1 and gID not in alreadyChecked:                         # an already indexed point belonging to a grain not yet tested?
          alreadyChecked[gID] = True                                        # remember not to check again
          disorientation = o.disorientation(orientations[gID])              # compare against that grain's orientation
          if disorientation.quaternion.w > cos_disorientation and \
             disorientation.quaternion.w >= bestDisorientation.quaternion.w:  # within disorientation threshold and better than current best?
            matched = True
            matchedID = gID                                                 # remember that grain
#            file['croak'].write('%i %f '%(matchedID,disorientation.quaternion.w))

            bestDisorientation = disorientation

    if not matched:                                                         # no match -> new grain found
      memberCounts += [1]                                                   # start new membership counter
      orientations += [o]                                                   # initialize with current orientation
      matchedID = g
      g += 1                                                                # increment grain counter
#      file['croak'].write('+')

    else:                                                                   # did match existing grain
      memberCounts[matchedID] += 1
#      file['croak'].write('got back %s is close by %f to %s\n'%(np.degrees(bestQ.asEulers()),np.degrees(2*np.arccos(bestDisorientation.quaternion.w)),np.degrees(bestFormerQ.asEulers())))
#      file['croak'].write('.%i %s'%(matchedID, orientations[matchedID-1].quaternion))
#      M = (1. - 1./memberCounts[matchedID-1]) * bestFormerQ.asM() + 1./memberCounts[matchedID-1] * bestQ.asM() # 4x4 matrix holding weighted quaternion outer products per grain
#      w,v = np.linalg.eigh(M)
#      avgQ = damask.Orientation(quaternion=v[:,w.argmax()],symmetry=options.symmetry)
#      file['croak'].write('new avg has misori of %f\n'%np.degrees(2*np.arccos(orientations[matchedID-1].disorientation(avgQ)[0].quaternion.w)))
#      orientations[matchedID-1].quaternion = damask.Quaternion(v[:,w.argmax()])
#      orientations[matchedID-1] = damask.Orientation(quaternion = bestDisorientation.quaternion**(1./memberCounts[matchedID-1]) \
#                                                                * orientations[matchedID-1].quaternion,
#                                                     symmetry = options.symmetry)   # adjust average orientation taking newest member into account
#      file['croak'].write(' stored --> %s\n'%(np.degrees(orientations[matchedID-1].quaternion.asEulers())))
#      file['croak'].write('.')

    grainID[p] = matchedID                                                  # remember grain index assigned to point
    p += 1                                                                  # increment point

  bg.set_message('identifying similar orientations among %i grains...'%(len(orientations)))

  memberCounts = np.array(memberCounts)
  similarOrientations = [[] for i in xrange(len(orientations))]

  for i,orientation in enumerate(orientations):                             # compare each identified orientation...
    for j in xrange(i+1,len(orientations)):                                 # ...against all others that were defined afterwards
      if orientation.disorientation(orientations[j]).quaternion.w > cos_disorientation:      # similar orientations in both grainIDs?
        similarOrientations[i].append(j)                                    # remember in upper triangle...
        similarOrientations[j].append(i)                                    # ...and lower triangle of matrix

    if similarOrientations[i] != []:
      bg.set_message('grainID %i is as: %s'%(i,' '.join(map(lambda x:str(x),similarOrientations[i]))))

  stillShifting = True
  while stillShifting:
    stillShifting = False
    tick = time.clock()

    for p,gID in enumerate(grainID):                                        # walk through all points
      if p > 0 and p % 1000 == 0:

        time_delta = (time.clock()-tick) * (len(grainID) - p) / p
        bg.set_message('(%02i:%02i:%02i) shifting ID of point %i out of %i (grain count %i)...'%(time_delta//3600,time_delta%3600//60,time_delta%60,p,len(grainID),len(orientations)))
      if similarOrientations[gID] != []:                                    # orientation of my grainID is similar to someone else?
        similarNeighbors = defaultdict(int)                                 # dict holding frequency of neighboring grainIDs that share my orientation (freq info not used...)
        for i in kdtree.query_ball_point(kdtree.data[p],options.radius):    # check all neighboring points
          if grainID[i] in similarOrientations[gID]:                        # neighboring point shares my orientation?
            similarNeighbors[grainID[i]] += 1                               # remember its grainID
        if similarNeighbors != {}:                                          # found similar orientation(s) in neighborhood
          candidates = np.array([gID]+similarNeighbors.keys())              # possible replacement grainIDs for me
          grainID[p] = candidates[np.argsort(memberCounts[candidates])[-1]] # adopt ID that is most frequent in overall dataset
          memberCounts[gID]        -= 1                                     # my former ID loses one fellow
          memberCounts[grainID[p]] += 1                                     # my new ID gains one fellow
          bg.set_message('%i:%i --> %i'%(p,gID,grainID[p]))                 # report switch of grainID
          stillShifting = True

  table.data_rewind()
  p = 0
  while table.data_read():                                                  # read next data line of ASCII table
    table.data_append(1+grainID[p])                                         # add grain ID
    table.data_write()                                                      # output processed line
    p += 1

  bg.set_message('done after %i seconds'%(time.clock()-start))

#  for i,o in enumerate(orientations):                                       # croak about average grain orientations
#    file['croak'].write('%i: %s\n'%(i,' '.join(map(str,o.quaternion.asEulers()))))

# ------------------------------------------ output result ---------------------------------------

  table.output_flush()                                                      # just in case of buffered ASCII table
  table.close()                                                             # close ASCII tables
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                             # overwrite old one with tmp new
