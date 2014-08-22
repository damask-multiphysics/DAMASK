#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,itertools,numpy,damask
from collections import defaultdict
from optparse import OptionParser, Option

scriptID = '$Id$'
scriptName = scriptID.split()[1]

# -----------------------------
class extendableOption(Option):
# -----------------------------
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



#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
synonyms = {
        'grid':   ['resolution'],
        'size':   ['dimension'],
          }
identifiers = {
        'grid':   ['a','b','c'],
        'size':   ['x','y','z'],
        'origin': ['x','y','z'],
          }
mappings = {
        'grid':            lambda x: int(x),
        'size':            lambda x: float(x),
        'origin':          lambda x: float(x),
        'homogenization':  lambda x: int(x),
        'microstructures': lambda x: int(x),
          }


parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Create seed file by taking microstructure indices from given ASCIItable column.
White and black-listing of microstructure indices is possible.

Examples:
--white 1,2,5 --index grainID isolates grainID entries of value 1, 2, and 5;
--black 1 --index grainID takes all grainID entries except for value 1.
""" + string.replace(scriptID,'\n','\\n')
)


parser.add_option('-p', '--positions',   dest = 'pos', type = 'string',
                                    help = 'coordinate label')
parser.add_option('--boundingbox',  dest = 'box', type = 'float', nargs = 6,
                                    help = 'min (x,y,z) and max (x,y,z) to specify bounding box [auto]')
parser.add_option('-i', '--index',  dest = 'index', type = 'string',
                                    help = 'microstructure index label')
parser.add_option('-w','--white',   dest = 'whitelist', action = 'extend', type = 'string', \
                                    help = 'white list of microstructure indices', metavar = '<LIST>')
parser.add_option('-b','--black',   dest = 'blacklist', action = 'extend', type = 'string', \
                                    help = 'black list of microstructure indices', metavar = '<LIST>')

parser.set_defaults(pos = 'pos')
parser.set_defaults(index = 'microstructure')
parser.set_defaults(box = [])
parser.set_defaults(whitelist = [])
parser.set_defaults(blacklist = [])

(options,filenames) = parser.parse_args()

datainfo = {                                                               # list of requested labels per datatype
             'scalar':     {'len':1,
                            'label':[]},
             'vector':     {'len':3,
                            'label':[]},
           }

if options.pos   != None:  datainfo['vector']['label'] += [options.pos]
if options.index != None:  datainfo['scalar']['label'] += [options.index]
options.whitelist = map(int,options.whitelist)
options.blacklist = map(int,options.blacklist)

#--- setup file handles --------------------------------------------------------------------------   
files = []
if filenames == []:
  files.append({'name':'STDIN',
                'input':sys.stdin,
                'output':sys.stdout,
                'croak':sys.stderr,
               })
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name,
                    'input':open(name),
                    'output':open(os.path.splitext(name)[0]+'.seeds','w'),
                    'croak':sys.stdout,
                    })

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  theTable = damask.ASCIItable(file['input'],file['output'],buffered = False)
  theTable.head_read()

# --------------- figure out columns to process
  active = defaultdict(list)
  column = defaultdict(dict)

  for datatype,info in datainfo.items():
    for label in info['label']:
      foundIt = False
      for key in ['1_'+label,label]:
        if key in theTable.labels:
          foundIt = True
          active[datatype].append(label)
          column[datatype][label] = theTable.labels.index(key)                 # remember columns of requested data
      if not foundIt:
        file['croak'].write('column %s not found...\n'%label)
        break


# ------------------------------------------ process data ---------------------------------------  

  theTable.data_readArray(list(itertools.chain.from_iterable(map(lambda x:[x+i for i in range(datainfo['vector']['len'])],
                                                                 [column['vector'][label] for label in active['vector']]))) + 
                       [column['scalar'][label] for label in active['scalar']])

#--- finding bounding box ------------------------------------------------------------------------------------
  boundingBox = numpy.array((numpy.amin(theTable.data[:,0:3],axis = 0),numpy.amax(theTable.data[:,0:3],axis = 0)))
  if len(options.box) == 6:
    boundingBox[0,:] = numpy.minimum(options.box[0:3],boundingBox[0,:])
    boundingBox[1,:] = numpy.maximum(options.box[3:6],boundingBox[1,:])

#--- rescaling coordinates ------------------------------------------------------------------------------------
  theTable.data[:,0:3] -= boundingBox[0,:]
  theTable.data[:,0:3] /= boundingBox[1,:]-boundingBox[0,:]


#--- filtering of grain voxels ------------------------------------------------------------------------------------
  mask = numpy.logical_and(\
         numpy.ones_like(theTable.data[:,3],bool) \
          if options.whitelist == [] \
          else              numpy.in1d(theTable.data[:,3].ravel(), options.whitelist).reshape(theTable.data[:,3].shape),
         numpy.ones_like(theTable.data[:,3],bool) \
          if options.blacklist == [] \
          else numpy.invert(numpy.in1d(theTable.data[:,3].ravel(), options.blacklist).reshape(theTable.data[:,3].shape))
          )
  theTable.data = theTable.data[mask]

# ------------------------------------------ output result ---------------------------------------  

# ------------------------------------------ assemble header ---------------------------------------  

  theTable.info = [
                   scriptID,
                   'size %s'%(' '.join(list(itertools.chain.from_iterable(zip(['x','y','z'],
                                                                          map(str,boundingBox[1,:]-boundingBox[0,:])))))),
                  ]
  theTable.labels_clear()
  theTable.labels_append(['x','y','z','microstructure'])                                  # implicitly switching label processing/writing on
  theTable.head_write()
  
  theTable.data_writeArray()
  theTable.output_flush()
  
  table.input_close()                                                       # close input ASCII table
  if file['name'] != 'STDIN':
    table.output_close()                                                    # close output ASCII table
