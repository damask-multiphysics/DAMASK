#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

# Makes postprocessing routines acessible from everywhere.
import os
from damask import Environment

damaskEnv = Environment()
baseDir = damaskEnv.relPath('processing/')
codeDir = damaskEnv.relPath('code/')
try:
  binDir = damaskEnv.options['DAMASK_BIN']
except:
  root=os.access('/usr/local/bin', os.W_OK)
  if root:
    binDir = '/usr/local/bin'
  else:
    binDir = os.path.join(os.getenv('HOME'),'bin')

if not os.path.isdir(binDir):
  os.mkdir(binDir)

#define ToDo list
bin_link = ['pre','post','misc']
            
for myDir in bin_link:
  myDir = os.path.abspath(os.path.join(baseDir,myDir))
  for myFile in os.listdir(myDir):
    src = os.path.abspath(os.path.join(myDir,myFile))
    sym_link = os.path.abspath(os.path.join(binDir,os.path.splitext(myFile)[0]))
    print sym_link,'-->',src
    if os.path.lexists(sym_link):
      os.remove(sym_link)    
    os.symlink(src,sym_link)            



