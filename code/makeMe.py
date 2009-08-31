#!/usr/bin/env python
# $Id$
import os,sys

sys.argv += ['' for i in range(2 - len(sys.argv))]
architectures = { 
                 'marc': { 
                          'parent': 'mpie_cpfem_marc.f90', 
                          'versions' : ['%%MARCVERSION%%','2007r1','2008r1'],
                          'substitutions' : {'%%REVISION%%': sys.argv[1],},
                         }, 
                }

for arch in architectures:
	me = architectures[arch]
	try:
		parentFile = open(me['parent'])
		parentContent = parentFile.readlines()
		parentFile.close()
	except IOError:
		print 'unable to open',me['parent']
		continue

	
	for version in me['versions'][1:]:
		childFile = open(version.join(os.path.splitext(me['parent'])),'w')
		for line in parentContent:
			for substitution in me['substitutions']:
				line = line.replace(substitution,me['substitutions'][substitution])
			childFile.write(line.replace(me['versions'][0],version))
		childFile.close()
