#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math,numpy
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP


# -----------------------------
class extendedOption(Option):
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


# ----------------------- MAIN -------------------------------

identifiers = {
        'resolution': ['a','b','c'],
        'dimension':  ['x','y','z'],
        'origin':     ['x','y','z'],
          }
mappings = {
        'resolution': lambda x: int(x),
        'dimension':  lambda x: float(x),
        'origin':     lambda x: float(x),
        'homogenization': lambda x: int(x),
        'maxgraincount': lambda x: int(x),
          }


parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
compress geometry files with ranges "a to b" and/or multiples "n of x".
""" + string.replace('$Id$','\n','\\n')
)

(options, filenames) = parser.parse_args()

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w')})

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': print file['name']

  #  get labels by either read the first row, or - if keyword header is present - the last line of the header

  firstline = file['input'].readline()
  m = re.search('(\d+)\s*head', firstline.lower())
  if m:
    headerlines = int(m.group(1))
    headers  = [firstline]+[file['input'].readline() for i in range(headerlines)]
  else:
    headerlines = 1
    headers = firstline

  content = file['input'].readlines()
  file['input'].close()

  info = {'resolution': [0,0,0],
          'dimension':  [0.0,0.0,0.0],
          'origin':     [0.0,0.0,0.0],
          'homogenization': 1,
          'maxgraincount': 0,
         }

  new_header = []
  for header in headers:
    headitems = map(str.lower,header.split())
    if headitems[0] in mappings.keys():
      if headitems[0] in identifiers.keys():
        for i in xrange(len(identifiers[headitems[0]])):
          info[headitems[0]][i] = \
            mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])
      else:
        info[headitems[0]] = mappings[headitems[0]](headitems[1])

  if info['resolution'] == [0,0,0]:
    print 'no resolution info found.'
    continue
  if info['dimension'] == [0.0,0.0,0.0]:
    print 'no dimension info found.'
    continue

  if file['name'] != 'STDIN':
    print 'resolution: %s'%(' x '.join(map(str,info['resolution'])))
    print 'dimension:  %s'%(' x '.join(map(str,info['dimension'])))
    print 'origin:     %s'%(' : '.join(map(str,info['origin'])))

  new_header.append("resolution\ta %i\tb %i\tc %i\n"%( 
    info['resolution'][0],
    info['resolution'][1],
    info['resolution'][2],))
  new_header.append("dimension\tx %f\ty %f\tz %f\n"%(
    info['dimension'][0],
    info['dimension'][1],
    info['dimension'][2]))
  new_header.append("origin\tx %f\ty %f\tz %f\n"%(
    info['origin'][0],
    info['origin'][1],
    info['origin'][2]))
  new_header.append("homogenization\t%i\n"%info['homogenization'])
  new_header.append("maxGrainCount\t%i\n"%info['maxgraincount'])

# ------------------------------------------ assemble header ---------------------------------------  

  file['output'].write('%i\theader\n'%(len(new_header))+''.join(new_header))

# ------------------------------------------ pack input ---------------------------------------  

  type = ''
  former = -1
  start = -1
  reps = 0

  for line in content:
    for current in map(int,line.split()):
      if current == former+1 and start+reps == former+1:   
        type = 'to'
        reps += 1
      elif current == former and start == former:
        type = 'of'
        reps += 1
      else:
        output = {'': '',
                  '.':  str(former)+'\n',
                  'to': '%i to %i\n'%(former-reps+1,former),
                  'of': '%i of %i\n'%(reps,former),
                  }[type]
        file['output'].write(output)
        type = '.'
        start = current
        reps = 1

      former = current
      
# write out last item...

  output = {'.':  str(former),
            'to': '%i to %i'%(former-reps+1,former),
            'of': '%i of %i'%(reps,former),
            }[type]
  file['output'].write(output+'\n')

# ------------------------------------------ output finalization ---------------------------------------  

  if file['name'] != 'STDIN':
    file['output'].close()
    os.rename(file['name']+'_tmp',file['name'])
