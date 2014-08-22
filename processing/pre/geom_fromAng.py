#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math,string,numpy as np
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP

scriptID = '$Id$'
scriptName = scriptID.split()[1]

#--------------------------------------------------------------------------------------------------
class extendableOption(Option):
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
parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Generate geometry description and material configuration from EBSD data in given square-gridded 'ang' file.
Two phases can be discriminated based on threshold value in a given data column.
""" + string.replace(scriptID,'\n','\\n')
)


parser.add_option('--column',              dest='column', type='int', metavar = 'int', \
                  help='data column to discriminate phase 1 from 2 [%default]')
parser.add_option('-t','--threshold',      dest='threshold', type='float', metavar = 'float', \
                  help='threshold value to discriminate phase 1 from 2 [%default]')
parser.add_option('--homogenization',      dest='homogenization', type='int', metavar = 'int', \
                  help='homogenization index to be used [%default]')
parser.add_option('--phase',               dest='phase', type='int', nargs = 2, metavar = 'int int', \
                  help='two phase indices to be used %default')
parser.add_option('--crystallite',         dest='crystallite', type='int', metavar = 'int', \
                  help='crystallite index to be used [%default]')
parser.add_option('-c', '--configuration', dest='config', action='store_true', \
                  help='output material configuration [%default]')
parser.add_option('--compress',            dest='compress', action='store_true', \
                  help='search for matching mircrostructure and texture and lump them [%default]')
parser.add_option('-a', '--axes',         dest='axes', type='string', nargs = 3, metavar = 'string string string', \
                  help='axes assignement of eulerangles x,y,z = %default')
    
                  
parser.set_defaults(column = 11)
parser.set_defaults(threshold = 0.5)
parser.set_defaults(homogenization = 1)
parser.set_defaults(phase          = [1,2])
parser.set_defaults(crystallite    = 1)
parser.set_defaults(config = False)
parser.set_defaults(compress= False)
parser.set_defaults(axes           = ['y','x','-z'])
(options,filenames) = parser.parse_args()

for i in options.axes:
  if i.lower() not in ['x','+x','-x','y','+y','-y','z','+z','-z']:
    file['croak'].write('invalid axes %s %s %s' %(options.axes[0],options.axes[1],options.axes[2]))
    sys.exit()

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
  if file['name'] != 'STDIN': file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')

  info = {
          'grid':   np.ones (3,'i'),
          'size':   np.zeros(3,'d'),
          'origin': np.zeros(3,'d'),
          'microstructures': 0,
          'homogenization':  options.homogenization
         }

  step           = [0,0]
  point          = 0
  for line in file['input']:
    words = line.split()
    if words[0] == '#':                                                                             # process initial comments block
      if len(words) > 2:
        if words[2].lower() == 'hexgrid': 
          file['croak'].write('The file has HexGrid format. Please first convert to SquareGrid...\n')
          sys.exit()
        if words[1] == 'XSTEP:':     step[0] = float(words[2])
        if words[1] == 'YSTEP:':     step[1] = float(words[2])
        if words[1] == 'NCOLS_ODD:':
          info['grid'][0] = int(words[2])
          eulerangles = np.zeros((info['grid'][0]*info['grid'][1],3),dtype='f')
          phase = np.zeros(info['grid'][0]*info['grid'][1],dtype='i')
        if words[1] == 'NROWS:':         
          info['grid'][1] = int(words[2])
          eulerangles = np.zeros((info['grid'][0]*info['grid'][1],3),dtype='f')
          phase = np.zeros(info['grid'][0]*info['grid'][1],dtype='i')
    else:                                                                                           # finished with comments block
      phase[point] = options.phase[{True:0,False:1}[float(words[options.column-1])<options.threshold]]
      eulerangles[point,...] = map(lambda x: float(x)*180.0/math.pi, words[:3])
      point += 1
 
  if info['grid'].prod() != point:
    file['croak'].write('Error: found %s microstructures. Header info in ang file might be wrong.\n'%point)
    sys.exit()

  if options.compress:
    texture = []
    microstructure = []
    otherPoint=-1                                                                                   # ensure to create first microstructure
    matPoints = np.zeros(info['grid'][0]*info['grid'][1],dtype='i')
    for myPoint in xrange(info['grid'][0]*info['grid'][1]):
      myTexture=-1
      for otherPoint in xrange(len(microstructure)):
        otherEulers = eulerangles[texture[microstructure[otherPoint][0]]]
        otherPhase  = phase[[microstructure[otherPoint][1]]]
        if all(eulerangles[myPoint]==otherEulers) and phase[myPoint] == otherPhase:  # common microstructure
           matPoints[myPoint] = otherPoint+1                                         # use other points microstructure
           otherPoint =-4                                                            # in no case, create new microstructure
           break
        elif all(eulerangles[myPoint]==otherEulers):                                 # found common texture and store it
           myTexture = microstructure[otherPoint][0]
      if otherPoint == len(microstructure)-1 or otherPoint == -2:                    #
        if myTexture == -1:
          texture.append(myPoint)
          myTexture = len(texture)-1
        microstructure.append([myTexture,myPoint])
        matPoints[myPoint] = len(microstructure)    # use the new microstructure

  else:
    texture = [i for i in xrange(info['grid'][0]*info['grid'][1])]
    microstructure = [[i+1,phase[i]] for i in xrange(info['grid'][0]*info['grid'][1])]
    matPoints = np.arange(info['grid'][0]*info['grid'][1],dtype='i')
    
  formatOut = 1+int(math.log10(len(texture)))
  textureOut =['\n\n<texture>']
  for i in xrange(len(texture)):
    textureOut +=       ['[Texture%s]\n'%str(i+1).zfill(formatOut) + \
                          'axes %s %s %s\n'%(options.axes[0],options.axes[1],options.axes[2]) +\
                          '(gauss)\tphi1 %4.2f\tPhi %4.2f\tphi2 %4.2f\tscatter 0.0\tfraction 1.0\n'%tuple(eulerangles[texture[i],...])
                         ]
  formatOut = 1+int(math.log10(len(microstructure)))
  microstructureOut =['<microstructure>']
  for i in xrange(len(microstructure)):
    microstructureOut += ['[Grain%s]\n'%str(i+1).zfill(formatOut) + \
                         'crystallite\t%i\n'%options.crystallite + \
                         '(constituent)\tphase %i\ttexture %s\tfraction 1.0\n'%(phase[microstructure[i][1]],microstructure[i][0]+1)
                         ]

  info['microstructures'] = len(microstructure)
  info['size'] = step[0]*info['grid'][0],step[1]*info['grid'][1],min(step)

#--- report ---------------------------------------------------------------------------------------
  file['croak'].write('grid     a b c:  %s\n'%(' x '.join(map(str,info['grid']))) + \
                      'size     x y z:  %s\n'%(' x '.join(map(str,info['size']))) + \
                      'origin   x y z:  %s\n'%(' : '.join(map(str,info['origin']))) + \
                      'homogenization:  %i\n'%info['homogenization'] + \
                      'microstructures: %i\n\n'%info['microstructures'])

  if np.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    sys.exit()
  if np.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    sys.exit()


#--- write data -----------------------------------------------------------------------------------
  if options.config:
    file['output'].write('\n'.join(microstructureOut) + \
                         '\n'.join(textureOut))
  else:
    header = [scriptID + ' ' + ' '.join(sys.argv[1:]) + '\n']
    header.append("grid\ta %i\tb %i\tc %i\n"%(info['grid'][0],info['grid'][1],info['grid'][2],))
    header.append("size\tx %f\ty %f\tz %f\n"%(info['size'][0],info['size'][1],info['size'][2],))
    header.append("origin\tx %f\ty %f\tz %f\n"%(info['origin'][0],info['origin'][1],info['origin'][2],))
    header.append("microstructures\t%i\n"%info['microstructures'])
    header.append("homogenization\t%i\n"%info['homogenization'])
    file['output'].write('%i\theader\n'%(len(header))+''.join(header))
    if options.compress:
      matPoints = matPoints.reshape((info['grid'][1],info['grid'][0]))
      np.savetxt(file['output'],matPoints,fmt='%0'+str(1+int(math.log10(np.amax(matPoints))))+'d')
    else:
      file['output'].write("1 to %i\n"%(info['microstructures']))
  
#--- output finalization -------------------------------------------------------------------------- 
  if file['name'] != 'STDIN':
    table.output_close()  
    os.rename(file['name']+'_tmp',os.path.splitext(file['name'])[0] + \
                                  {True: '_material.config',
                                   False:'.geom'}[options.config])
