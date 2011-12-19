#!/usr/bin/env python
# $Id$
# Writes version specific files for different MARC releases
import os,sys,string,re,damask

architectures = { 
                 'marc': { 
                          'parent': 'DAMASK_marc.f90', 
                          'versions' : ['%%MARCVERSION%%','2007r1','2008r1','2010','2011'],
                         }, 
                }

bin_link = { \
        '' : [
                'DAMASK_spectral.exe',
             ],
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

makefile = open(os.path.join(baseDir,'makefile'))
content = makefile.readlines()
makefile.close()
makefile = open(os.path.join(baseDir,'makefile'),'w')
for line in content:
  m = re.match(r'(FFTW|ACML)ROOT\s*:?=',line)
  if m: line = '%sROOT := %s\n'%(m.group(1),damaskEnv.pathInfo[m.group(1).lower()])
  makefile.writelines(line)
makefile.close()

# compiling spectral code
compile = raw_input("Do you want to compile the spectral code now? (y/n) ")
if (compile == 'y' or compile == 'Y'):
  compiler_switches = raw_input("Please give compiling switches (Enter to use default) ")
  os.system('make --directory %s clean'%(baseDir))
  os.system('make --directory %s %s'%(baseDir,compiler_switches))

if '--clean' in [s.lower() for s in sys.argv]:
  os.system('make --directory %s clean'%baseDir)


for dir in bin_link:
  for file in bin_link[dir]:
    src = os.path.abspath(os.path.join(baseDir,dir,file))
    if (file == ''):
      sym_link = os.path.abspath(os.path.join(damaskEnv.binDir(),dir))
    else:
      sym_link = os.path.abspath(os.path.join(damaskEnv.binDir(),os.path.splitext(file)[0]))
    print sym_link,'-->',src
    if os.path.lexists(sym_link):
      os.remove(sym_link)    
    os.symlink(src,sym_link)
