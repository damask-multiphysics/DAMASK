#!/usr/bin/env python
# $Id$
# Writes version specific files for different MARC releases
import os,sys,string,damask_tools

architectures = { 
                 'marc': { 
                          'parent': 'DAMASK_marc.f90', 
                          'versions' : ['%%MARCVERSION%%','2007r1','2008r1','2010'],
                         }, 
                }

bin_link = { \
        '' : [
                'DAMASK_spectral.exe',
             ],
            }
            
damask_variables = damask_tools.DAMASK_TOOLS()
baseDir = damask_variables.relPath('code/')

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

# changing dirs in make file

makefile = open(os.path.join(baseDir,'makefile'))
content = makefile.readlines()
makefile.close()
makefile = open(os.path.join(baseDir,'makefile'),'w')
for line in content:
  if line.startswith('FFTWPATH'):
    line='FFTWPATH =%s\n'%(damask_variables.pathInfo['fftw'])
    print line
  if line.startswith('ACMLROOT'):
    line='ACMLROOT =%s\n'%(damask_variables.pathInfo['acml'])
    print line
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
      sym_link = os.path.abspath(os.path.join(damask_variables.binDir(),dir))
    else:
      sym_link = os.path.abspath(os.path.join(damask_variables.binDir(),os.path.splitext(file)[0]))
    print sym_link,'-->',src
    if os.path.lexists(sym_link):
      os.remove(sym_link)    
    os.symlink(src,sym_link)
