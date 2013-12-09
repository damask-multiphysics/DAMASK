#!/usr/bin/env python

import os,sys,string,re,math,numpy
import damask
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP

scriptID = '$Id$'
scriptName = scriptID.split()[1]

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


parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
Add RGB color value corresponding to TSL-OIM scheme for inverse pole figures.
""" + string.replace(scriptID,'\n','\\n')
)

parser.add_option('-p', '--pole', dest='pole', type='float', nargs=3,\
                  help = 'lab frame direction for inverse pole figure %default')
parser.add_option('-s', '--symmetry', dest='symmetry', type='string', \
                  help = 'crystal symmetry [%default]')
parser.add_option('-e', '--eulers',   dest='eulers', type='string', \
                  help = 'Euler angles label')
parser.add_option('-m', '--matrix',   dest='matrix', type='string', \
                  help = 'orientation matrix label')
parser.add_option('-a',               dest='a', type='string', \
                  help = 'crystal frame a vector label')
parser.add_option('-b',               dest='b', type='string', \
                  help = 'crystal frame b vector label')
parser.add_option('-c',               dest='c', type='string', \
                  help = 'crystal frame c vector label')
parser.add_option('-q', '--quaternion', dest='quaternion', type='string', \
                  help = 'quaternion label')

parser.set_defaults(pole = [0.0,0.0,1.0])
parser.set_defaults(symmetry = 'cubic')

(options, filenames) = parser.parse_args()

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

pole = numpy.array(options.pole)
pole /= numpy.linalg.norm(pole)

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w'), 'croak':sys.stderr})

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  table = damask.ASCIItable(file['input'],file['output'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info
  table.info_append(string.replace(scriptID,'\n','\\n') + '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out columns to process
  active = {}
  column = {}
  head = []

  for datatype,info in datainfo.items():
    for label in info['label']:
      foundIt = False
      for key in ['1_'+label,label]:
        if key in table.labels:
          foundIt = True
          if datatype not in active: active[datatype] = []
          if datatype not in column: column[datatype] = {}
          active[datatype].append(label)
          column[datatype][label] = table.labels.index(key)                 # remember columns of requested data
      if not foundIt:
        file['croak'].write('column %s not found...\n'%label)
        break

  table.labels_append(['%i_IPF_%g%g%g'%(i+1,options.pole[0],options.pole[1],options.pole[2]) for i in xrange(3)])

# ------------------------------------------ assemble header ---------------------------------------  

  table.head_write()

# ------------------------------------------ process data ---------------------------------------  

  while table.data_read():                                                  # read next data line of ASCII table

    if input == 'eulers':
      o = damask.Orientation(Eulers=numpy.array(map(float,table.data[column['vector'][options.eulers]:\
                                                                      column['vector'][options.eulers]+datainfo['vector']['len']])),
                             symmetry=options.symmetry).reduced()
    elif input == 'matrix':
      o = damask.Orientation(matrix=numpy.array([map(float,table.data[column['tensor'][options.matrix]:\
                                                                      column['tensor'][options.matrix]+datainfo['tensor']['len']])]),
                             symmetry=options.symmetry).reduced()
    elif input == 'frame':
      o = damask.Orientation(matrix=numpy.array([map(float,table.data[column['vector'][options.a]:\
                                                                      column['vector'][options.a]+datainfo['vector']['len']] + \
                                                           table.data[column['vector'][options.b]:\
                                                                      column['vector'][options.b]+datainfo['vector']['len']] + \
                                                           table.data[column['vector'][options.c]:\
                                                                      column['vector'][options.c]+datainfo['vector']['len']]
                                                    )]).reshape(3,3),
                             symmetry=options.symmetry).reduced()
    elif input == 'quaternion':
      o = damask.Orientation(quaternion=numpy.array(map(float,table.data[column['quaternion'][options.quaternion]:\
                                                                          column['quaternion'][options.quaternion]+datainfo['quaternion']['len']])),
                             symmetry=options.symmetry).reduced()

    table.data_append(o.IPFcolor(pole))

    table.data_write()                                                      # output processed line

# ------------------------------------------ output result ---------------------------------------  

  table.output_flush()                                                      # just in case of buffered ASCII table

  if file['name'] != 'STDIN':
    file['input'].close()                                                   # close input ASCII table
    file['output'].close()                                                  # close output ASCII table
    os.rename(file['name']+'_tmp',file['name'])                             # overwrite old one with tmp new
