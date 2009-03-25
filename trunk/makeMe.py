#!/usr/bin/env python

import os,sys

architectures = { 
                 'marc': { 
                          'parent': 'mpie_cpfem_marc.f90', 
                          'versions' : ['%%MARCVERSION%%','2007r1','2008r1'], 
                         }, 
                }

for arch in architectures:
	try:
		parent = architectures[arch]['parent']
		parentFile = open(parent)
		parentContent = parentFile.readlines()
		parentFile.close()
	except IOError:
		print 'unable to open',parent
		continue

	for version in architectures[arch]['versions'][1:]:
		childFile = open(os.path.splitext(parent)[0]+version+os.path.splitext(parent)[1],'w')
		for line in parentContent:
			childFile.write(line.replace(architectures[arch]['versions'][0],version))
		childFile.close()
