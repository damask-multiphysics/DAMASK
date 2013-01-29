#!/usr/bin/env python
import os,string,re,damask

bin_link = { \
            '.' : [
                    'DAMASK_spectral.exe',
                  ],
           }

damaskEnv = damask.Environment('../../')          # script location relative to root
baseDir = damaskEnv.relPath('code/')

for dir in bin_link:
  for file in bin_link[dir]:
    print baseDir
    print dir
    print file
    src = os.path.abspath(os.path.join(baseDir,dir,file))
    print src
    os.path.exists(src)
    if os.path.exists(src): 
      sym_link = os.path.abspath(os.path.join(damaskEnv.binDir(),\
                                              {True: dir,
                                               False:os.path.splitext(file)[0]}[file == '']))
      if os.path.lexists(sym_link):
        os.remove(sym_link)
      os.symlink(src,sym_link)
      print sym_link,'-->',src
