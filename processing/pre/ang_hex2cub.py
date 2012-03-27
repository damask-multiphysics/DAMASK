#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import string,os,sys
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

parser = OptionParser(option_class=extendedOption, usage='%prog [geomfile[s]]', description = """
Converts EBSD Data stored in *.ang files from hex to cub

""" + string.replace('$Id: spectral_geomCheck.py 1084 2011-11-09 15:37:45Z MPIE\c.zambaldi $','\n','\\n')
)

(options, filenames) = parser.parse_args()

# ------------------------------------------ setup file handles ---------------------------------------  

columnX   = 3 # 0,1,2,3 (python notation!)
counterX  = 0
counterY  = 0
addPoints = 0

files = []
for name in filenames:
  if os.path.exists(name):
    files.append(   {'name':name, 'input':open(name),'output':open('cub_'+name, 'w')})

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  print file['name']
  for line in file['input']:
    if line.split()[0]=='#':
      if len(line.split())>2:                                                                               # possibly interesting information
        if line.split()[2]=='SqrGrid': 
          print 'The file is already a square grid file.'
          sys.exit()
        if line.split()[2]=='HexGrid': line='# GRID: SqrGrid\n'
        if line.split()[1]=='XSTEP:':      stepSizeX = float(line.split()[2]) 
        if line.split()[1]=='YSTEP:':      stepSizeY = float(line.split()[2]) 
        if line.split()[1]=='NCOLS_EVEN:': NColsEven = int(line.split()[2]) 
      file['output'].write(line)
    else:                                                                                                   # finished reading of header
      lineSplit=line.split()
      x = float(lineSplit[columnX])
      y = float(lineSplit[columnX+1])
      lineFirstPart =''
      lineLastPart =''
      for i in xrange(columnX):
        lineFirstPart =lineFirstPart+' '+lineSplit[i]
      for i in xrange(len(lineSplit)- (columnX+2)):
        lineLastPart =lineLastPart+' '+lineSplit[i+columnX+2]

      file['output'].write(lineFirstPart+' '+\
                   str((counterX+addPoints)*stepSizeY)+' '+str(y)+' '+\
                   lineLastPart+'\n')
      if x + stepSizeX - (counterX+addPoints+1)*stepSizeY > 0.5*stepSizeY:                                 # double point (interpolation error)
        addPoints+=1
        file['output'].write(lineFirstPart+' '+\
                   str((counterX+addPoints)*stepSizeY)+' '+str(y)+' '+\
                   lineLastPart+'\n')

      if(counterX == NColsEven + counterY%2):                                                              # new row (odd and even differ by 1)                   
        counterY+=1
        counterX=0
        addPoints=0
      counterX+=1
