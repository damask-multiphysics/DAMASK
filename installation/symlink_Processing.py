#!/usr/bin/env python2
# -*- coding: UTF-8 no BOM -*-

# Makes postprocessing routines acessible from everywhere.
import os,sys
import damask

damaskEnv = damask.Environment()
baseDir = damaskEnv.relPath('processing/')
binDir = damaskEnv.options['DAMASK_BIN']

if not os.path.isdir(binDir):
  os.mkdir(binDir)

#define ToDo list
processing_subDirs = ['pre','post','misc',]
processing_extensions = ['.py','.sh',]
            
for subDir in processing_subDirs:
  theDir = os.path.abspath(os.path.join(baseDir,subDir))

  for theFile in os.listdir(theDir):
    if os.path.splitext(theFile)[1] in processing_extensions:                           # only consider files with proper extensions

      src      = os.path.abspath(os.path.join(theDir,theFile))
      sym_link = os.path.abspath(os.path.join(binDir,os.path.splitext(theFile)[0]))

      if os.path.lexists(sym_link):
        os.remove(sym_link)
        sys.stdout.write(sym_link)
      else:
        sys.stdout.write(damask.util.emph(sym_link))

      os.symlink(src,sym_link)
      sys.stdout.write(' -> '+src+'\n')
