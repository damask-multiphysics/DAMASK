#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os,sys,math,string,re,numpy, damask 
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
  

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
identifiers = {
        'resolution': ['a','b','c'],
        'dimension':  ['x','y','z'],
        'origin':     ['x','y','z'],
          }
mappings = {
        'grains': lambda x: int(x),
        'resolution': lambda x: int(x),
        'origin': lambda x: float(x),
        'dimension': lambda x: float(x),
        'homogenization': lambda x: int(x),
          }

parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
Generate geometry description and material configuration by standard Voronoi tessellation of given seeds file.
""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                  help='output geom file with two-dimensional data arrangement')

parser.set_defaults(twoD   = False)

(options,filenames) = parser.parse_args()


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

  firstline = file['input'].readline()
  m = re.search('(\d+)\s*head', firstline.lower())
  if m:
    headerlines = int(m.group(1))
    headers  = [firstline]+[file['input'].readline() for i in range(headerlines)]
  else:
    headerlines = 1
    headers = firstline


  info = {'grains': 0,
          'resolution': numpy.array([0,0,0]),
          'dimension': numpy.array([1.0,1.0,1.0]),
          'origin':     numpy.array([0.0,0.0,0.0]),
          'homogenization': 1,
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
  
  content = []
  for line in file['input']:
    for word in  line.split(): content.append(word)
  file['input'].close()
  content=numpy.reshape(numpy.tile(numpy.reshape(content,info['resolution']),(1,1,1,3,3,3)),
               info['resolution'].prod()*27)

  info['grains'] =  max(numpy.array(content,'i'))
  file['croak'].write('grains:         %i\n'%info['grains'] + \
                      'resolution:     %s\n'%(' x '.join(map(str,info['resolution']))) + \
                      'dimension:      %s\n'%(' x '.join(map(str,info['dimension']))) + \
                      'origin:         %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization: %i\n'%info['homogenization'])

  new_header.append("grains\t%i\n"%( info['grains']))
  new_header.append("resolution\ta %i\tb %i\tc %i\n"%( 
    info['resolution'][0]*3,
    info['resolution'][1]*3,
    info['resolution'][2]*3,))
  new_header.append("dimension\tx %f\ty %f\tz %f\n"%(
    info['dimension'][0],
    info['dimension'][1],
    info['dimension'][2],))
  new_header.append("origin\tx %f\ty %f\tz %f\n"%(
    info['origin'][0],
    info['origin'][1],
    info['origin'][2],))
  new_header.append("homogenization\t%i\n"%info['homogenization'])

# -------------------------------------- write data to file  ----------------------------------

  formatwidth = 1+int(math.log10(info['grains']))
  file['output'].write('%i\theader\n'%(len(new_header)) + ''.join(new_header))
  
  for n in xrange(info['resolution'][1:3].prod()*9):   
    file['output'].write({ True: ' ',
                           False:'\n'}[options.twoD].\
                           join(map(lambda x: str(x).rjust(formatwidth),\
                           content[n*3*info['resolution'][0]:(n+1)*3*info['resolution'][0]]))+'\n')
  
# ------------------------------------------ output finalization ---------------------------------------  

  if file['name'] != 'STDIN':
    file['output'].close()
    os.rename(file['name']+'_tmp',os.path.splitext(file['name'])[0]+'_27periodicCopies.geom')
