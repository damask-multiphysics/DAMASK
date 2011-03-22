#!/usr/bin/env python
# $Id$
# Writes version specific files for different MARC releases
import os,sys

architectures = { 
                 'marc': { 
                          'parent': 'mpie_cpfem_marc.f90', 
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
