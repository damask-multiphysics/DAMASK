#!/usr/bin/env python

import os,sys,math,string,numpy
from optparse import OptionParser, Option

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


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Generate geometry description and material configuration from EBSD data in given square-gridded 'ang' file.
Two phases can be discriminated based on threshold value in a given data column.
""" + string.replace('$Id$','\n','\\n')
)


parser.add_option('--column',          dest='column', type='int', \
                  help='data column to discriminate phase 1 from 2 [%default]')
parser.add_option('-t','--threshold',  dest='threshold', type='float', \
                  help='threshold value to discriminate phase 1 from 2 [%default]')
parser.add_option('--homogenization', dest='homogenization', type='int', \
                  help='homogenization index to be used [%default]')
parser.add_option('--phase', dest='phase', type='int', nargs = 2, \
                  help='phase indices to be used %default')
parser.add_option('--crystallite', dest='crystallite', type='int', \
                  help='crystallite index to be used [%default]')
parser.add_option('-c', '--configuration', dest='config', action='store_true', \
                  help='output material configuration [%default]')
                  
parser.set_defaults(column = 11)
parser.set_defaults(threshold = 0.5)
parser.set_defaults(homogenization = 1)
parser.set_defaults(phase          = [1,2])
parser.set_defaults(crystallite    = 1)
parser.set_defaults(config = False)

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

  point          = 0
  step           = [0,0]
  grid     = [1,1]
  microstructure = ['<microstructure>']
  texture        = ['<texture>']

  for line in file['input']:
    words = line.split()
    if words[0] == '#':                                                 # process initial comments block
      if len(words) > 2:
        if words[1] == 'HexGrid': 
          file['croak'].write('The file has HexGrid format. Please first convert to SquareGrid...\n'); break
        if words[1] == 'XSTEP:':     step[0] = float(words[2])
        if words[1] == 'YSTEP:':     step[1] = float(words[2])
        if words[1] == 'NCOLS_ODD:': 
          grid[0] =   int(words[2]); formatwidth = 1+int(math.log10(grid[0]*grid[1]))
        if words[1] == 'NROWS:':     
          grid[1] =   int(words[2]); formatwidth = 1+int(math.log10(grid[0]*grid[1]))
    else:                                                               # finished with comments block
      if options.config:                                                # write configuration (line by line)
        point += 1
        me = str(point).zfill(formatwidth)
        microstructure += ['[Grain%s]\n'%me + \
                           'crystallite\t%i\n'%options.crystallite + \
                           '(constituent)\tphase %i\ttexture %s\tfraction 1.0\n'%(options.phase[{True:0,False:1}[float(words[options.column-1])<options.threshold]],me)
                          ]
        texture +=        ['[Grain%s]\n'%me + \
                           '(gauss)\tphi1 %4.2f\tPhi %4.2f\tphi2 %4.2f\tscatter 0.0\tfraction 1.0\n'%tuple(map(lambda x: float(x)*180.0/math.pi, words[:3]))
                          ]
      else:
        grid.append(1)
        file['output'].write("6 header\n" + \
                             "$Id$ \n"
                             "grid\ta %i\tb %i\tc 1\n"%(grid[0],grid[1]) + \
                             "size\tx %g\ty %g\tz %g\n"%(step[0]*grid[0],step[1]*grid[1],min(step)) + \
                             "origin\tx 0\ty 0\tz 0\n" + \
                             "microstructures\t%i\n"%(grid[0]*grid[1]) + \
                             "homogenization %i\n"%options.homogenization + \
                             "1 to %i\n"%(grid[0]*grid[1]))
        break
  file['croak'].write('grid     a b c:  %s\n'%(' x '.join(map(str,grid))) + \
                      'size     x y z:  %s\n'%(' x '.join(map(str,[step[0]*grid[0],step[1]*grid[1],min(step)]))) + \
                      'origin   x y z:  %s\n'%(' : '.join(map(str,[0.0,0.0,0.0]))) + \
                      'microstructures: %i\n'%(grid[0]*grid[1]) + \
                      'homogenization:  %i\n'%options.homogenization)
  
  if options.config:
    file['output'].write('\n'.join(microstructure) + \
                         '\n'.join(texture))
  
  # ------------------------------------------ output finalization ---------------------------------------  

  if file['name'] != 'STDIN':
    file['output'].close()
    os.rename(file['name']+'_tmp',os.path.splitext(file['name'])[0] + \
                                  {True: '_material.config',
                                   False:'.geom'}[options.config])
