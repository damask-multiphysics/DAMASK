#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

# Makes postprocessing routines acessible from everywhere.
import os,sys
from damask import Environment

BOLD = '\033[1m'
ENDC = '\033[0m'

damaskEnv = Environment()
baseDir = damaskEnv.relPath('processing/')
codeDir = damaskEnv.relPath('code/')
try:
  binDir = damaskEnv.options['DAMASK_BIN']
except:
  binDir = '/usr/local/bin' if os.access('/usr/local/bin', os.W_OK) else os.path.join(os.getenv('HOME'),'bin')

if not os.path.isdir(binDir):
  os.mkdir(binDir)

#define ToDo list
processing_subDirs = ['pre','post','misc',]
processing_extensions = ['.py',]
            
for subDir in processing_subDirs:
  theDir = os.path.abspath(os.path.join(baseDir,subDir))

  for theFile in os.listdir(theDir):
    if os.path.splitext(theFile)[1] in processing_extensions:                           # omit anything not fitting our script extensions (skip .py.bak, .py~, and the like)

      src      = os.path.abspath(os.path.join(theDir,theFile))
      sym_link = os.path.abspath(os.path.join(binDir,os.path.splitext(theFile)[0]))

      if os.path.lexists(sym_link):
        os.remove(sym_link)
        sys.stdout.write(sym_link)
      else:
        sys.stdout.write(BOLD + sym_link + ENDC)

      os.symlink(src,sym_link)
      sys.stdout.write(' -> '+src+'\n')
