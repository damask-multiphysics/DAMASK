#!/usr/bin/env python
import os,site,sys

# More elegant was to determine DAMASK install dir?
#p=os.path.abspath(__file__)
#l=p.split('/')
#p=p[:-len(l[-1])]
#print p
#site.addsitedir(p) # adds the paths in the *.pth file to the python search path

basepath=os.path.expanduser('~')+'/DAMASK/'


sys.path.insert(0,basepath+'processing')
sys.path.insert(0,basepath+'processing/post')
sys.path.insert(0,basepath+'processing/pre')
sys.path.insert(0,basepath+'processing/setup')
#print sys.path



