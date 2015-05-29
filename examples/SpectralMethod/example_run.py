#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-
import sys

resolutions = [16,32,64]
resolution = resolutions[0]

try:
  resolution = int(sys.argv[1])
except:
  pass

if resolution not in resolutions:
  resolution = resolutions[0]
  
from subprocess import call
call('make run%s'%('x'.join([str(resolution)]*3)), shell=True)