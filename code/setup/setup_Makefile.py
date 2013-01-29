#!/usr/bin/env python
import os,string,re,damask


architectures = { 
                 'marc': { 
                          'parent': 'DAMASK_marc.f90', 
                          'versions' : ['%%MARCVERSION%%','2010','2011','2012'],
                         }, 
                }


damaskEnv = damask.Environment('../../')          # script location relative to root
baseDir = damaskEnv.relPath('code/')

for arch in architectures:
  me = architectures[arch]
  try:
    parentFile = open(baseDir+os.sep+me['parent'])
    parentContent = parentFile.readlines()
    parentFile.close()
  except IOError:
    print 'unable to open',me['parent']
    continue

  
  for version in me['versions'][1:]:
    childFile = open(baseDir+os.sep+version.join(os.path.splitext(me['parent'])),'w')
    for line in parentContent:
      childFile.write(line.replace(me['versions'][0],version))
    childFile.close()

# changing dirs in makefile
makefile = open(os.path.join(baseDir,'Makefile'))
content = makefile.readlines()
makefile.close()
makefile = open(os.path.join(baseDir,'Makefile'),'w')
for line in content:
  m = re.match(r'(FFTW|IMKL|ACML|LAPACK)ROOT\s*\?=',line)
  if m:
    if m.group(1).lower() in damaskEnv.pathInfo:
      substitution = damaskEnv.pathInfo[m.group(1).lower()]
    else:
      substitution = ''
    line = '%sROOT ?= %s\n'%(m.group(1),substitution)
  makefile.write(line)
makefile.close()
















