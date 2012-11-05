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
          }
mappings = {
        'resolution': lambda x: int(x),
        'grains':     lambda x: int(x),
          }

parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
generates geom file and material_config file using seeds file

""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-r', '--resolution', dest='resolution', type='int', nargs = 3, \
                                       help='a,b,c resolution of specimen')
parser.add_option('-d', '--dimension', dest='dimension', type='float', nargs = 3, \
                                       help='x,y,z dimension of specimen')
parser.add_option('--homogenization', dest='homogenization', type='int', \
                                      help='homogenization index to be used')
parser.add_option('--phase', dest='phase', type='int', \
                             help='phase index to be used')
parser.add_option('-c', '--configuration', dest='config', action='store_true', \
                                           help='output material configuration')
parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                  help='output geom file with two-dimensional data arrangement')

                                   
parser.set_defaults(resolution = [0,0,0])
parser.set_defaults(dimension  = [0.0,0.0,0.0])
parser.set_defaults(homogenization = 1)
parser.set_defaults(phase          = 1)
parser.set_defaults(config = False)
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

  content = file['input'].readlines()
  file['input'].close()

  info = {'grains': 0,
          'resolution': numpy.array([0,0,0]),
          'dimension':  numpy.array(options.dimension),
          'origin':     numpy.array([0.0,0.0,0.0]),
          'homogenization': options.homogenization,
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

  if info['grains'] == 0:
    file['croak'].write('no grains found.\n')
    continue
  if info['grains'] != len(content):
    file['croak'].write('grain data not matching grain count...\n')
    info['grains'] = min(info['grains'],len(content))

  if 0 not in options.resolution:                                    # user-specified resolution
    info['resolution'] = numpy.array(options.resolution)

  if info['resolution'].all() == 0:
    file['croak'].write('no resolution info found.\n')
    continue

  twoD = info['resolution'][2] < 2

  for i in xrange(3):
    if info['dimension'][i] <= 0.0:                                 # any invalid dimension?
      info['dimension'][i] = float(info['resolution'][i])/max(info['resolution'])
      file['croak'].write('rescaling dimension %i...\n'%i)

  file['croak'].write('grains:         %i\n'%info['grains'] + \
                      'resolution:     %s\n'%(' x '.join(map(str,info['resolution']))) + \
                      'dimension:      %s\n'%(' x '.join(map(str,info['dimension']))) + \
                      'origin:         %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization: %i\n'%info['homogenization'])

  new_header.append("resolution\ta %i\tb %i\tc %i\n"%( 
    info['resolution'][0],
    info['resolution'][1],
    info['resolution'][2],))
  new_header.append("dimension\tx %f\ty %f\tz %f\n"%(
    info['dimension'][0],
    info['dimension'][1],
    info['dimension'][2],))
  new_header.append("origin\tx %f\ty %f\tz %f\n"%(
    info['origin'][0],
    info['origin'][1],
    info['origin'][2],))
  new_header.append("homogenization\t%i\n"%info['homogenization'])

# -------------------------------------- prepare data ----------------------------------

  formatwidth = 1+int(math.log10(info['grains']))
  coords = numpy.zeros((3,info['grains']),'d')
  eulers = numpy.zeros((3,info['grains']),'d')

  for i in xrange(info['grains']):
    coords[:,i] = map(float,content[i].split()[:3])*info['dimension']
    eulers[:,i] = map(float,content[i].split()[3:6])
 
# -------------------------------------- switch according to task ----------------------------------

  if options.config:
    file['output'].write('<microstructure>\n')
    for i in xrange(info['grains']):
      file['output'].write('\n[Grain%s]\n'%(str(i+1).zfill(formatwidth)) + \
                           'crystallite 1\n' + \
                           '(constituent)\tphase %i\ttexture %s\tfraction 1.0\n'%(options.phase,str(i+1).rjust(formatwidth)))
  
    file['output'].write('\n<texture>\n')
    for i in xrange(info['grains']):
      file['output'].write('\n[Grain%s]\n'%(str(i+1).zfill(formatwidth)) + \
                           '(gauss)\tphi1 %g\tPhi %g\tphi2 %g\tscatter 0.0\tfraction 1.0\n'%(eulers[0,i],eulers[1,i],eulers[2,i]))
    
  else:
    file['output'].write('%i\theader\n'%(len(new_header)) + ''.join(new_header))
  
    N = info['resolution'].prod()
    shift = 0.5*info['dimension']/info['resolution']                            # shift by half of side length to center of element
    undeformed = numpy.zeros((3,N),'d')

    for i in xrange(N):
      undeformed[0,i] = info['dimension'][0]\
                       * float(i                                              % info['resolution'][0])\
                                                                         /float(info['resolution'][0])
      undeformed[1,i] = info['dimension'][1]\
                      * float(i//info['resolution'][0]                        % info['resolution'][1])\
                                                                         /float(info['resolution'][1])
      undeformed[2,i] = info['dimension'][2]\
                      * float(i//info['resolution'][0]//info['resolution'][1] % info['resolution'][2])\
                                                                         /float(info['resolution'][2])
      undeformed[:,i] += shift
      
    indices = damask.core.math.math_nearestNeighborSearch(3,\
              numpy.array(([1.0,0.0,0.0],\
                           [0.0,1.0,0.0],\
                           [0.0,0.0,1.0]),'d'),\
              info['dimension'],\
              N,info['grains'],undeformed,coords)//3**3 + 1                     # floor division to kill periodic images
  
    for n in xrange(info['resolution'][1:3].prod()):                            # loop over 2nd and 3rd dimension
      file['output'].write({ True: ' ',
                             False:'\n'}[options.twoD].\
                             join(map(lambda x: str(x).rjust(formatwidth),\
                                      indices[n*info['resolution'][0]:(n+1)*info['resolution'][0]]))+'\n')
  
    missing = 0
    for i in xrange(info['grains']):
      if i+1 not in indices: missing += 1
    file['croak'].write({True:'all',False:'only'}[missing == 0] + ' %i grains mapped.\n'%(info['grains']-missing))
  
# ------------------------------------------ output finalization ---------------------------------------  

  if file['name'] != 'STDIN':
    file['output'].close()
    os.rename(file['name']+'_tmp',os.path.splitext(file['name'])[0] + \
                                  {True: '_material.config',
                                   False:'.geom'}[options.config])

