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
        'maxmicrostructure': lambda x: int(x),
          }


parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
Unpack geometry files containing ranges "a to b" and/or "n of x" multiples (exclusively in one line).
""" + string.replace('$Id: spectral_geomCanvas.py 1576 2012-06-26 18:08:50Z MPIE\p.eisenlohr $','\n','\\n')
)

parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                  help='output geom file with two-dimensional data arrangement')

parser.set_defaults(twoD = False)

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
          'maxmicrostructure': 0,
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

  format = {True:  info['resolution'][0],
            False: 1}[options.twoD]

  if file['name'] != 'STDIN':
    print 'resolution:    %s'%(' x '.join(map(str,info['resolution'])))
    print 'dimension:     %s'%(' x '.join(map(str,info['dimension'])))
    print 'origin:        %s'%(' : '.join(map(str,info['origin'])))

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
  new_header.append("maxMicrostructure\t%i\n"%info['maxmicrostructure'])
  if info['maxmicrostructure'] != 0:
    digits = 1+int(math.log10(int(info['maxmicrostructure'])))
  else:
    digits = 1+int(math.log10(int(info['resolution'][0]*info['resolution'][1]*info['resolution'][2])))
  print digits
# ------------------------------------------ assemble header ---------------------------------------  

  file['output'].write('%i\theader\n'%(len(new_header))+''.join(new_header))

# ------------------------------------------ unpack input ---------------------------------------  

  wordsWritten = 0
  for line in content:
    words = map(str.lower,line.split())
    if len(words) > 1:        # any packing keywords?
      if (words[1] == 'to'): words = map(str,range(int(words[0]),int(words[2])+1))
      if (words[1] == 'of'): words = [words[2]]*int(words[0])

    for word in words:
      wordsWritten += 1
      file['output'].write(word.zfill(digits)+{True:'\n',False:' '}[wordsWritten%format == 0])    # newline every format words

# ------------------------------------------ output finalization ---------------------------------------  

  if file['name'] != 'STDIN':
    file['output'].close()
    os.rename(file['name']+'_tmp',file['name'])
