#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,sys
import damask

MarcReleases =[ \
               '2014',
               '2014.2',
               '2015',
               '2016'
              ]

damaskEnv = damask.Environment()
baseDir = damaskEnv.relPath('code/')
binDir  = damaskEnv.options['DAMASK_BIN']

if not os.path.isdir(binDir):
  os.mkdir(binDir)

sys.stdout.write('\nMSC.Marc versioning...\n\n')
theMaster = 'DAMASK_marc.f90'

for version in MarcReleases:
  src = os.path.abspath(os.path.join(baseDir,theMaster))
  if os.path.exists(src): 
    sym_link = os.path.abspath(os.path.join(baseDir,'DAMASK_marc{}.f90'.format(version)))                    
    if os.path.lexists(sym_link):
      os.remove(sym_link)
      output = version
    else:
      output = damask.util.emph(version)

    sys.stdout.write(' '+output+'\n')
    os.symlink(theMaster,sym_link)
