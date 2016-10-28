#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

import os,sys
import damask

bin_link = { \
            '.' : [
                    'DAMASK_spectral.exe',
                    'DAMASK_FEM.exe'
                  ],
           }

MarcReleases =[ \
               '2011',
               '2012',
               '2013',
               '2013.1',
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

sys.stdout.write('\nsymbolic linking...\n')

for subDir in bin_link:
  theDir = os.path.abspath(os.path.join(baseDir,subDir))
  sys.stdout.write('\n'+binDir+' ->\n'+theDir+damask.util.deemph(' ...')+'\n')

  for theFile in bin_link[subDir]:
    theName,theExt = os.path.splitext(theFile)
    src = os.path.abspath(os.path.join(theDir,theFile))

    if os.path.exists(src): 
      sym_link = os.path.abspath(os.path.join(binDir,subDir if theFile == '' else theName))

      if os.path.lexists(sym_link):
        os.remove(sym_link)
        output = theName+damask.util.deemph(theExt)
      else:
        output = damask.util.emph(theName)+damask.util.deemph(theExt)

      sys.stdout.write(damask.util.deemph('... ')+output+'\n')
      os.symlink(src,sym_link)


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
