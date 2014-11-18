#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import string,os,sys
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Converts ang files (EBSD Data) from hexagonal grid to a pixel grid

""", version = scriptID)

parser.add_option('-x', dest='columnX', type='int', metavar='int', \
                        help='column containing x coordinates [%default]')

parser.set_defaults(columnX = 3)

(options,filenames) = parser.parse_args()

counterX  = 0
counterY  = 0
addPoints = -1                                                                                      # No of doubled points (must be the same for each odd/even line, initializing with -1 make countin easy!)

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
                    'output':open(os.path.splitext(name)[0]+'_cub'+os.path.splitext(name)[1], 'w'),
                    'croak':sys.stdout,
                 })

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  file['croak'].write('\033[1m' + scriptName + '\033[0m: ' + (file['name'] if file['name'] != 'STDIN' else '') + '\n')

  x = 0
  for line in file['input']:
    lineSplit=line.split()
    
    if lineSplit[0]=='#':
      if len(lineSplit)>2:                                                                          # possibly interesting information
        if line.split()[2]=='SqrGrid': 
          file['croak'].write('The file is already a square grid file.')
          sys.exit()
        if lineSplit[1]=='XSTEP:':      stepSizeX = float(lineSplit[2]) 
        if lineSplit[1]=='YSTEP:':      stepSizeY = float(lineSplit[2])
        
        if lineSplit[2]=='HexGrid':     line='# GRID: SqrGrid\n'                                    # comments are not read by OIM, but still better to be correct
        if lineSplit[1]=='NCOLS_ODD:':
          NCols = int(int(lineSplit[2])*stepSizeX/stepSizeY)
          line='# NCOLS_ODD: %d\n'% NCols
        if lineSplit[1]=='NCOLS_EVEN:':
          line='# NCOLS_EVEN: %d\n'% NCols
          
      file['output'].write(line)
    else:                                                                                           # finished reading of header      
      xOld = x
      x = float(lineSplit[options.columnX])                                                         # current (original) x positions
      
      if x >  xOld:                                                                                 # same line, increase X
        counterX+=1
      else:                                                                                         # new line, increase in Y, reset X
        counterY+=1
        addPoints = -1                                                                              # to start at zero
        counterX=0
      
      lineFirstPart =''                                                                             # split line around x and y coordinate
      for i in xrange(options.columnX):
        lineFirstPart =lineFirstPart+' '+lineSplit[i]
      lineLastPart =''
      for i in xrange(len(lineSplit)- (options.columnX+2)):
        lineLastPart =lineLastPart+' '+lineSplit[i+options.columnX+2]
        
      if counterX+addPoints < NCols:
        file['output'].write('%s %.6f %.6f  %s\n' %(lineFirstPart,(counterX+addPoints)*stepSizeY,   # write with new x and y position 
                                                                  counterY*stepSizeY,lineLastPart))                     

      if x - (counterX+addPoints)*stepSizeY > 0.5*stepSizeY and counterX+addPoints+1 < NCols:       # double point (interpolation error)
                   
        addPoints+=1
        file['output'].write('%s %.6f %.6f  %s\n' %(lineFirstPart,(counterX+addPoints)*stepSizeY,\
                                                                   counterY*stepSizeY,lineLastPart))
