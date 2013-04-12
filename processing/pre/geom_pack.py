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
        'grid':    ['a','b','c'],
        'size':    ['x','y','z'],
        'origin':  ['x','y','z'],
          }
mappings = {
        'grid':            lambda x: int(x),
        'size':            lambda x: float(x),
        'origin':          lambda x: float(x),
        'homogenization':  lambda x: int(x),
        'microstructures': lambda x: int(x),
          }


parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
compress geometry files with ranges "a to b" and/or multiples "n of x".
""" + string.replace('$Id$','\n','\\n')
)

(options, filenames) = parser.parse_args()

# ------------------------------------------ setup file handles ---------------------------------------  

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
                    'output':open(name+'_tmp','w'),
                    'croak':sys.stdout,
                    })

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': file['croak'].write(file['name']+'\n')

#  get labels by either read the first row, or - if keyword header is present - the last line of the header
  firstline = file['input'].readline()
  m = re.search('(\d+)\s*head', firstline.lower())
  if m:
    headerlines = int(m.group(1))
    headers  = [file['input'].readline() for i in range(headerlines)]
  else:
    headerlines = 1
    headers = firstline

  content = file['input'].readlines()
  file['input'].close()

  info = {'grid':           [0,0,0],
          'size':           [0.0,0.0,0.0],
          'origin':         [0.0,0.0,0.0],
          'homogenization':  1,
          'microstructures': 0,
         }

  new_header = []
  for header in headers:
    headitems = map(str.lower,header.split())
    if headitems[0] == 'resolution': headitems[0] = 'grid'
    if headitems[0] == 'dimension':  headitems[0] = 'size'
    if headitems[0] in mappings.keys():
      if headitems[0] in identifiers.keys():
        for i in xrange(len(identifiers[headitems[0]])):
          info[headitems[0]][i] = \
            mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])
      else:
        info[headitems[0]] = mappings[headitems[0]](headitems[1])
    new_header.append(header)

  if info['grid'] == [0,0,0]:
    file['croak'].write('no grid info found.\n')
    continue
  if info['size'] == [0.0,0.0,0.0]:
    file['croak'].write('no size info found.\n')
    continue

  file['croak'].write('grid     a b c:  %s\n'%(' x '.join(map(str,info['grid']))) + \
                      'size     x y z:  %s\n'%(' x '.join(map(str,info['size']))) + \
                      'origin   x y z:  %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization:  %i\n'%info['homogenization'] + \
                      'microstructures: %i\n'%info['microstructures'])

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
