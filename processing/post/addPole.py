#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add x,y coordinates of stereographic projection of given direction (pole) in crystal frame.

""", version = scriptID)

parser.add_option('-p', '--pole',       dest='pole', type='float', nargs=3, metavar='float float float',
                                        help = 'crystal frame direction for pole figure %default')
parser.add_option('--polar',            dest='polar', action='store_true',
                                        help = 'output polar coordinates r,phi [%default]')
parser.add_option('-e', '--eulers',     dest='eulers', metavar='string',
                                        help = 'Euler angles label')
parser.add_option('-d', '--degrees',    dest='degrees', action='store_true',
                                        help = 'Euler angles are given in degrees [%default]')
parser.add_option('-m', '--matrix',     dest='matrix', metavar='string',
                                        help = 'orientation matrix label')
parser.add_option('-a',                 dest='a', metavar='string',
                                        help = 'crystal frame a vector label')
parser.add_option('-b',                 dest='b', metavar='string',
                                        help = 'crystal frame b vector label')
parser.add_option('-c',                 dest='c', metavar='string',
                                        help = 'crystal frame c vector label')
parser.add_option('-q', '--quaternion', dest='quaternion', metavar='string',
                                        help = 'quaternion label')
parser.set_defaults(pole = [1.0,0.0,0.0])
parser.set_defaults(degrees = False)
parser.set_defaults(polar   = False)

(options, filenames) = parser.parse_args()

datainfo = {                                                                                        # list of requested labels per datatype
             'tensor':     {'len':9,
                            'label':[]},
             'vector':     {'len':3,
                            'label':[]},
             'quaternion': {'len':4,
                            'label':[]},
           }

input = []
if options.eulers     != None:
  datainfo['vector']['label'] += [options.eulers]
  input += ['eulers']
if options.a          != None and \
   options.b          != None and \
   options.c          != None:
  datainfo['vector']['label'] += [options.a,options.b,options.c]
  input += ['frame']
if options.matrix     != None:
  datainfo['tensor']['label'] += [options.matrix]
  input += ['matrix']
if options.quaternion != None:
  datainfo['quaternion']['label'] += [options.quaternion]
  input += ['quaternion']

if len(input) != 1: parser.error('needs exactly one input format...')
input = input[0]

toRadians = np.pi/180.0 if options.degrees else 1.0                                               # rescale degrees to radians
pole = np.array(options.pole)
pole /= np.linalg.norm(pole)

# --- loop over input files -------------------------------------------------------------------------
if filenames == []:
  filenames = ['STDIN']

for name in filenames:
  if name == 'STDIN':
    file = {'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m\n')
  else:
    if not os.path.exists(name): continue
    file = {'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],buffered = False)                          # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

  column = {}
  missingColumns = False

  for datatype,info in datainfo.items():
    for label in info['label']:
      key = list(set([label, '1_'+label]) & set(table.labels))
      if key == []:
        file['croak'].write('column %s not found...\n'%label)
        missingColumns = True                                                                       # break if label not found
      else:
        column[label] = table.labels.index(key[0])                                                  # remember columns of requested data

  if missingColumns:
    continue

# ------------------------------------------ assemble header ---------------------------------------
  table.labels_append(['%i_pole_%g%g%g'%(i+1,options.pole[0],options.pole[1],options.pole[2]) for i in xrange(2)])
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    if input == 'eulers':
      o = damask.Orientation(Eulers=toRadians*\
                 np.array(map(float,table.data[column[options.eulers]:\
                                               column[options.eulers]+datainfo['vector']['len']])),
                             )
    elif input == 'matrix':
      o = damask.Orientation(matrix=\
                 np.array([map(float,table.data[column[options.matrix]:\
                                                column[options.matrix]+datainfo['tensor']['len']])]).reshape(np.sqrt(datainfo['tensor']['len']),
                                                                                                             np.sqrt(datainfo['tensor']['len'])).transpose(),
                             )
    elif input == 'frame':
      o = damask.Orientation(matrix=\
                 np.array([map(float,table.data[column[options.a]:\
                                                column[options.a]+datainfo['vector']['len']] + \
                                     table.data[column[options.b]:\
                                                column[options.b]+datainfo['vector']['len']] + \
                                     table.data[column[options.c]:\
                                                column[options.c]+datainfo['vector']['len']]
                                                    )]).reshape(3,3),
                             )
    elif input == 'quaternion':
      o = damask.Orientation(quaternion=\
                 np.array(map(float,table.data[column[options.quaternion]:\
                                               column[options.quaternion]+datainfo['quaternion']['len']])),
                             )

    rotatedPole = o.quaternion*pole                                                                 # rotate pole according to crystal orientation
    (x,y) = rotatedPole[0:2]/(1.+abs(pole[2]))                                                      # stereographic projection

    table.data_append([np.sqrt(x*x+y*y),np.arctan2(y,x)] if options.polar else [x,y])                                                                      # cartesian coordinates

    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table (works for stdin)
  table.output_close()                                                                              # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
