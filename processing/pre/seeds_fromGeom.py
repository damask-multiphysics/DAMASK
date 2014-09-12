#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,numpy,damask
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
Create seed file taking microstructure indices from given geom file but excluding black-listed grains.

Examples:
""" + string.replace(scriptID,'\n','\\n')
)


parser.add_option('-w','--white',   dest='whitelist', action='extend', \
                                    help='white list of grain IDs', metavar='<LIST>')
parser.add_option('-b','--black',   dest='blacklist', action='extend', \
                                    help='black list of grain IDs', metavar='<LIST>')

parser.set_defaults(whitelist = [])
parser.set_defaults(blacklist = [])

(options,filenames) = parser.parse_args()

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

  theTable = damask.ASCIItable(file['input'],file['output'],labels = False,buffered = False)
  theTable.head_read()

#--- interpret header ----------------------------------------------------------------------------
  info = {
          'grid':    numpy.zeros(3,'i'),
          'size':    numpy.zeros(3,'d'),
          'origin':  numpy.zeros(3,'d'),
          'homogenization':  0,
          'microstructures': 0,
         }
  newInfo = {
          'grid':    numpy.zeros(3,'i'),
          'origin':  numpy.zeros(3,'d'),
          'microstructures': 0,
         }
  extra_header = []

  for header in theTable.info:
    headitems = map(str.lower,header.split())
    if len(headitems) == 0: continue                                                              # skip blank lines
    for synonym,alternatives in synonyms.iteritems():
      if headitems[0] in alternatives: headitems[0] = synonym
    if headitems[0] in mappings.keys():
      if headitems[0] in identifiers.keys():
        for i in xrange(len(identifiers[headitems[0]])):
          info[headitems[0]][i] = \
            mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])
      else:
        info[headitems[0]] = mappings[headitems[0]](headitems[1])
    else:
      extra_header.append(header)

  file['croak'].write('grid     a b c:  %s\n'%(' x '.join(map(str,info['grid']))) + \
                      'size     x y z:  %s\n'%(' x '.join(map(str,info['size']))) + \
                      'origin   x y z:  %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization:  %i\n'%info['homogenization'] + \
                      'microstructures: %i\n'%info['microstructures'])

  if numpy.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    continue
  if numpy.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    continue
  if 'origin' not in info:
    info['origin'] = numpy.zeros(3)

#--- read data ------------------------------------------------------------------------------------
  microstructure = numpy.zeros(info['grid'].prod(),'i')                                            # initialize as flat array
  i = 0
  while theTable.data_read():
    items = theTable.data
    if len(items) > 2:
      if   items[1].lower() == 'of': items = [int(items[2])]*int(items[0])
      elif items[1].lower() == 'to': items = xrange(int(items[0]),1+int(items[2]))
      else:                            items = map(int,items)
    else:                              items = map(int,items)

    s = len(items)
    microstructure[i:i+s] = items
    i += s


# ------------------------------------------ assemble header ---------------------------------------  

  theTable.info = [
                   scriptID,
                   "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
                   "size\tx %i\ty %i\tz %i"%(info['size'][0],info['size'][1],info['size'][2],),
                   "origin\tx %i\ty %i\tz %i"%(info['origin'][0],info['origin'][1],info['origin'][2],),
                  ]
  theTable.labels_clear()
  theTable.labels_append(['x','y','z','microstructure'])                                  # implicitly switching label processing/writing on
  theTable.head_write()
  
#--- filtering of grain voxels ------------------------------------------------------------------------------------
  theTable.data_clear()
  i = 0
  outputDead = False
  coord = numpy.zeros(3,'d')
  for coord[2] in xrange(info['grid'][2]):
    for coord[1] in xrange(info['grid'][1]):
      for coord[0] in xrange(info['grid'][0]):
        if (options.whitelist == [] and options.blacklist == []) or \
           (options.whitelist != [] and microstructure[i]     in options.whitelist) or \
           (options.blacklist != [] and microstructure[i] not in options.blacklist):
          theTable.data = list((coord+0.5)/info['grid'])+[microstructure[i]]
          outputDead = not theTable.data_write()
        i += 1
        if outputDead: break
      if outputDead: break
    if outputDead: break

# ------------------------------------------ output result ---------------------------------------  

  outputDead or theTable.output_flush()                                        # just in case of buffered ASCII table

  table.input_close()                                                       # close input ASCII table
  if file['name'] != 'STDIN':
    table.output_close()                                                    # close output ASCII table