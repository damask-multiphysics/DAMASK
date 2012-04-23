#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,re,sys,math,string
from optparse import OptionParser, Option
    
def writeKey(outType, begin, counter):
  if outType == 'range':
    a=str(begin) +' to '+  str(begin+counter-1)+'\n'
    file['output'].write(a)
  if outType == 'copy':
    b=str(counter)+ ' copies of '+str(begin)+'\n' 
    file['output'].write(b)   
  if outType == 'none': 
    c=str(begin)+'\n'
    file['output'].write(c)
    
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
packing/unpacking geometry file 
""" + string.replace('$Id$','\n','\\n')
)


parser.add_option('-p','--pack',       dest='pack', action='store_true', \
                                       help='pack input file [%default]')
parser.add_option('-u','--unpack',     dest='unpack', action='store_true', \
                                       help='unpack input file [%default]')
parser.add_option('-l','--linelength', dest='lineLength',type='int', nargs=1, \
                                       help='length of line [%default]')
                                   

parser.set_defaults(pack = False)
parser.set_defaults(unpack = False)
parser.set_defaults(lineLength = 1)

(options,filenames) = parser.parse_args()

if options.pack and options.unpack:
  print 'cannot pack and unpack at the same time'
# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w+')})

# ------------------------------------------ loop over input files ---------------------------------------  
for file in files:
  counter = 1
  begin = -1
  outType = ''                                                                                        
  current=-1                                                                                        # current material ID
  charactersWritten=0                                                                               # number of words per line in the output file
  currentLine=0                   
  if file['name'] != 'STDIN': print file['name']
  for line in file['input']:
    words=line.split()
    currentLine+=1
    if('head' in line): headNumber=int(words[0])+1   
    if(currentLine<=headNumber): file['output'].write(line)
    else:
      if options.unpack:                                                                            # unpacking the geometry file
        words = line.split()
        if (words[1].lower()=='to'): line = ''.join(['%d '%i for i in xrange(int(words[0]), int(words[2])+1)])  
        if (words[1].lower()=='copies' and words[2].lower()=='of'): line=''.join(['%s '%words[3] for i in range(int(words[0]))])           
        words = line.split() 
        for i in xrange(len(words)):
          file['output'].write(words[i]+' ')
          charactersWritten+=1
          if(charactersWritten==options.lineLength):        
            file['output'].write('\n')
            charactersWritten=0
            
      if options.pack:                                                                              # packing the geometry file
        numbers=line.split()
        i=len(numbers)
        for x in xrange(i):
          last=current
          current = int(numbers[x])
          if current == last+1 and current==begin+counter:                                  
            counter+=1
            outType = 'range'
          elif current == begin:
            counter+=1
            outType = 'copy'
          else:
            writeKey(outType, begin, counter)
            counter = 1
            outType= 'none'
            begin = current
   
  writeKey(outType, begin, counter)                                                                 # just for packing (last line), outType not defined for unpacking
  os.rename(file['name']+'_tmp',file['name']) 
