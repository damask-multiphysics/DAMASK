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

wd = os.path.join(os.path.dirname(sys.argv[0]),'..')

for arch in architectures:
	me = architectures[arch]
	try:
		parentFile = open(wd+os.sep+me['parent'])
		parentContent = parentFile.readlines()
		parentFile.close()
	except IOError:
		print 'unable to open',me['parent']
		continue

	
	for version in me['versions'][1:]:
		childFile = open(wd+os.sep+version.join(os.path.splitext(me['parent'])),'w')
		for line in parentContent:
			childFile.write(line.replace(me['versions'][0],version))
		childFile.close()

# changing dirs in make file
damask_variables = damask_tools.DAMASK_TOOLS()
makefile = open(os.path.join(damask_variables.rootDir(),'code/makefile'))
content = makefile.readlines()
makefile.close()
makefile = open(os.path.join(damask_variables.rootDir(),'code/makefile'),'w')
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
if raw_input("Do you want to compile the spectral code now? (y/n) ") is 'y' or 'Y':
  compiler_switches = raw_input("Please give compiling switches (Enter to use default) ")
  os.system('make --directory %s clean'%(wd))
  print compiler_switches
  os.system('make --directory %s %s'%(wd,compiler_switches))

if '--clean' in [s.lower() for s in sys.argv]:
  os.system('make --directory %s clean'%wd)
