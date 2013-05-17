#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,sys,math,string,re,numpy, damask 
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP 

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
  

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
identifiers = {
        'grid':  ['a','b','c'],
          }
mappings = {
        'grid':    lambda x: int(x),
        'grains':  lambda x: int(x),
          }

parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
Generate geometry description and material configuration by standard Voronoi tessellation of given seeds file.
""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-g', '--grid', dest='grid', type='int', nargs = 3, \
                  help='a,b,c grid of hexahedral box [from seed file]')
parser.add_option('-s', '--size', dest='size', type='float', nargs = 3, \
                  help='x,y,z size of hexahedral box [1.0 along largest grid point number]')
parser.add_option('--homogenization', dest='homogenization', type='int', \
                  help='homogenization index to be used [%default]')
parser.add_option('--phase', dest='phase', type='int', \
                  help='phase index to be used [%default]')
parser.add_option('--crystallite', dest='crystallite', type='int', \
                  help='crystallite index to be used [%default]')
parser.add_option('-c', '--configuration', dest='config', action='store_true', \
                  help='output material configuration [%default]')
parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                  help='output geom file with two-dimensional data arrangement [%default]')

                                   
parser.set_defaults(grid = [0,0,0])
parser.set_defaults(size  = [0.0,0.0,0.0])
parser.set_defaults(homogenization = 1)
parser.set_defaults(phase          = 1)
parser.set_defaults(crystallite    = 1)
parser.set_defaults(config = False)
parser.set_defaults(twoD   = False)
                                   
(options,filenames) = parser.parse_args()

#--- setup file handles ---------------------------------------------------------------------------  
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


#--- loop over input files ------------------------------------------------------------------------
for file in files:
  if file['name'] != 'STDIN': file['croak'].write(file['name']+'\n')

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

#--- interprete header ----------------------------------------------------------------------------
  info = {
          'grid':    numpy.zeros(3,'i'),
          'size':    numpy.array(options.size),
          'origin':  numpy.zeros(3,'d'),
          'grains':  0,          
          'homogenization':  0,
         }

  new_header = []
  for header in headers:
    headitems = map(str.lower,header.split())
    if headitems[0] == 'resolution': headitems[0] = 'grid'
    if headitems[0] in mappings.keys():
      if headitems[0] in identifiers.keys():
        for i in xrange(len(identifiers[headitems[0]])):
          info[headitems[0]][i] = \
            mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])
      else:
        info[headitems[0]] = mappings[headitems[0]](headitems[1])

  if info['grains'] != len(content):
    file['croak'].write('grain data not matching grain count...\n')
    info['grains'] = min(info['grains'],len(content))
  
  if 0 not in options.grid:                                                                         # user-specified grid
    info['grid'] = numpy.array(options.grid)

  for i in xrange(3):
    if info['size'][i] <= 0.0:                                                                      # any invalid size?
      info['size'][i] = float(info['grid'][i])/max(info['grid'])
      file['croak'].write('rescaling size %i...\n'%i)

  file['croak'].write('grains to map:  %i\n'%info['grains'] + \
                      'grid     a b c: %s\n'%(' x '.join(map(str,info['grid']))) + \
                      'size     x y z: %s\n'%(' x '.join(map(str,info['size']))) + \
                      'origin   x y z: %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization: %i\n'%info['homogenization'])
  
  if numpy.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    sys.exit()
  if numpy.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    sys.exit()
  if info['grains'] == 0:
    file['croak'].write('no grain info found.\n')
    sys.exit()

#--- prepare data ---------------------------------------------------------------------------------
  formatwidth = 1+int(math.log10(info['grains']))
  coords = numpy.zeros((3,info['grains']),'d')
  eulers = numpy.zeros((3,info['grains']),'d')

  for i in xrange(info['grains']):
    coords[:,i] = map(float,content[i].split()[:3])*info['size']
    eulers[:,i] = map(float,content[i].split()[3:6])

#--- switch according to task ---------------------------------------------------------------------
  if options.config:                                                                                # write config file
    file['output'].write('<microstructure>\n')
    for i in xrange(info['grains']):
      file['output'].write('\n[Grain%s]\n'%(str(i+1).zfill(formatwidth)) + \
                           'crystallite %i\n'%options.crystallite + \
                           '(constituent)\tphase %i\ttexture %s\tfraction 1.0\n'%(options.phase,str(i+1).rjust(formatwidth)))
  
    file['output'].write('\n<texture>\n')
    for i in xrange(info['grains']):
      file['output'].write('\n[Grain%s]\n'%(str(i+1).zfill(formatwidth)) + \
                           '(gauss)\tphi1 %g\tPhi %g\tphi2 %g\tscatter 0.0\tfraction 1.0\n'%(eulers[0,i],eulers[1,i],eulers[2,i]))

  else:                                                                                             # write geometry file  
    twoD = info['grid'][2] < 2
    N = info['grid'].prod()
    shift = 0.5*info['size']/info['grid']                                                           # shift by half of side length to center of element
    undeformed = numpy.zeros((3,N),'d')

    for i in xrange(N):
      undeformed[0,i] = info['size'][0]\
                       * float(i                                              % info['grid'][0])\
                                                                         /float(info['grid'][0])
      undeformed[1,i] = info['size'][1]\
                      * float(i//info['grid'][0]                        % info['grid'][1])\
                                                                         /float(info['grid'][1])
      undeformed[2,i] = info['size'][2]\
                      * float(i//info['grid'][0]//info['grid'][1] % info['grid'][2])\
                                                                         /float(info['grid'][2])
      undeformed[:,i] += shift
      
    indices = damask.core.math.periodicNearestNeighbor(\
              info['size'],\
              numpy.eye(3),\
              undeformed,coords)//3**3 + 1                                                          # floor division to kill periodic images
    missing = 0
    for i in xrange(info['grains']):
      if i+1 not in indices: missing += 1
    file['croak'].write({True:'all',False:'only'}[missing == 0] + ' %i grains mapped.\n'%(info['grains']-missing))

#--- write header ---------------------------------------------------------------------------------
    new_header.append("$Id$ \n")
    new_header.append("grid\ta %i\tb %i\tc %i\n"%(info['grid'][0],info['grid'][1],info['grid'][2],))
    new_header.append("size\tx %f\ty %f\tz %f\n"%(info['size'][0],info['size'][1],info['size'][2],))
    new_header.append("origin\tx %f\ty %f\tz %f\n"%(info['origin'][0],info['origin'][1],info['origin'][2],))
    new_header.append("microstructures\t%i\n"%(info['grains']-missing))
    new_header.append("homogenization\t%i\n"%info['homogenization'])
    file['output'].write('%i\theader\n'%(len(new_header)) + ''.join(new_header))

# --- write microstructure information ------------------------------------------------------------
    for n in xrange(info['grid'][1:3].prod()):                                                        # loop over 2nd and 3rd size
      file['output'].write({ True: ' ',
                             False:'\n'}[options.twoD].\
                             join(map(lambda x: str(x).rjust(formatwidth),\
                                      indices[n*info['grid'][0]:(n+1)*info['grid'][0]]))+'\n')
   
#--- output finalization --------------------------------------------------------------------------
  if file['name'] != 'STDIN':
    file['output'].close()
    os.rename(file['name']+'_tmp',os.path.splitext(file['name'])[0] + \
                                  {True: '_material.config',
                                   False:'.geom'}[options.config])
