#!/usr/bin/env python
import os,sys
#import site
# maybe there is an alternative by using site.addsitedir() and *.pth file(s)?
damask_root=os.getenv('DAMASK_ROOT')
if damask_root is None:
  print('Environment variable DAMASK_ROOT not set.\nPlease run DAMASK/damask_env.sh first.')
  sys.exit()
print('Setting PYTHONPATH for DAMASK')  
path_list=['processing'
           ,'processing/setup'
           ,'processing/pre'
           ,'processing/post'
           ]
for p in path_list:
  if p not in sys.path:
    sys.path.insert(0, damask_root+'/'+p)
print sys.path



