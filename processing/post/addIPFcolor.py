#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add RGB color value corresponding to TSL-OIM scheme for inverse pole figures.

""", version = scriptID)

parser.add_option('-p', '--pole',       dest='pole', type='float', nargs=3, metavar='float float float',
                  help = 'lab frame direction for inverse pole figure [%default]')
parser.add_option('-s', '--symmetry',   dest='symmetry', type='choice', 
                  choices=damask.Symmetry.lattices[1:], metavar='string',
                  help = 'crystal symmetry [cubic] {%s} '%(', '.join(damask.Symmetry.lattices[1:])))
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

parser.set_defaults(pole = (0.0,0.0,1.0),
                    symmetry = 'cubic',
                    degrees = False,
                   )

(options, filenames) = parser.parse_args()

input = [options.eulers     != None,
         options.a          != None and \
         options.b          != None and \
         options.c          != None,
         options.matrix     != None,
         options.quaternion != None,
        ]

if np.sum(input) != 1: parser.error('needs exactly one input format...')

(label,dim,inputtype) = [(options.eulers,3,'eulers'),
                         ([options.a,options.b,options.c],[3,3,3],'frame'),
                         (options.matrix,9,'matrix'),
                         (options.quaternion,4,'quaternion'),
                        ][np.where(input)[0][0]]                                                                              # select input label that was requested
toRadians = math.pi/180.0 if options.degrees else 1.0                                               # rescale degrees to radians
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

  table = damask.ASCIItable(file['input'],file['output'],buffered=False)                            # make unbuffered ASCII_table
  table.head_read()                                                                                 # read ASCII header info

  if not np.all(table.label_dimension(label) == dim):
    file['croak'].write('input %s has wrong dimension %i...\n'%(label,dim))
    table.close(dismiss = True)                                                                     # close ASCIItable and remove empty file
    continue

  column = table.label_index(label)

# ------------------------------------------ assemble header ---------------------------------------
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append(['%i_IPF_%g%g%g_%s'%(i+1,options.pole[0],options.pole[1],options.pole[2],options.symmetry.lower()) for i in xrange(3)])
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    if inputtype == 'eulers':
      o = damask.Orientation(Eulers=toRadians*\
                             np.array(map(float,table.data[column:column+3])),
                             symmetry=options.symmetry).reduced()
    elif inputtype == 'matrix':
      o = damask.Orientation(matrix=\
                             np.array([map(float,table.data[column:column+9])]).reshape(3,3).transpose(),
                             symmetry=options.symmetry).reduced()
    elif inputtype == 'frame':
      o = damask.Orientation(matrix=\
                             np.array([map(float,table.data[column[0]:column[0]+3] + \
                                                 table.data[column[1]:column[1]+3] + \
                                                 table.data[column[2]:column[2]+3]
                                     )]).reshape(3,3),
                             symmetry=options.symmetry).reduced()
    elif inputtype == 'quaternion':
      o = damask.Orientation(quaternion=\
                             np.array(map(float,table.data[column:column+4])),
                             symmetry=options.symmetry).reduced()

    table.data_append(o.IPFcolor(pole))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output result -----------------------------------------
  outputAlive and table.output_flush()                                                              # just in case of buffered ASCII table

  table.input_close()                                                                               # close input ASCII table (works for stdin)
  table.output_close()                                                                              # close output ASCII table (works for stdout)
  if file['name'] != 'STDIN':
    os.rename(file['name']+'_tmp',file['name'])                                                     # overwrite old one with tmp new
