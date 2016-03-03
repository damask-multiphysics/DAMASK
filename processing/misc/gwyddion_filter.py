#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,string,scipy
import numpy as np
import damask
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Apply filter(s) to Gwyddion data.
""" + string.replace(scriptID,'\n','\\n')
)

for option in ['opening',
               'closing',
               'erosion',
               'dilation',
               'average',
               'median',
               ]:
  parser.add_option('-%s'%option[0], '--%s'%option,   dest=option, type='int',
                  help = 'stencil size for %s filter'%option)
  parser.set_default(option, 0)

(options, filenames) = parser.parse_args()


# ------------------------------------------ read Gwyddion data ---------------------------------------  

for file in filenames:
  filters = ''
  header = []
  with open(file,'r') as f:
    for line in f:
      pieces = line.split()
      if pieces[0] != '#': break
      if pieces[1] == 'Width:':  width  = float(pieces[2])
      if pieces[1] == 'Height:': height = float(pieces[2])
      header.append(line.lstrip('#').strip())
      
    elevation = np.loadtxt(file)#*1e6

    if options.opening > 0:
      elevation = scipy.ndimage.morphology.grey_opening(elevation,options.opening)
      filters += '_opening%i'%options.opening
    if options.closing > 0:
      elevation = scipy.ndimage.morphology.grey_closing(elevation,options.closing)
      filters += '_closing%i'%options.closing
    if options.erosion > 0:
      elevation = scipy.ndimage.morphology.grey_erosion(elevation,options.erosion)
      filters += '_erosion%i'%options.erosion
    if options.dilation > 0:
      elevation = scipy.ndimage.morphology.grey_dilation(elevation,options.dilation)
      filters += '_dilation%i'%options.dilation
    if options.average > 0:
      elevation = scipy.ndimage.filters.uniform_filter(elevation,options.average)
      filters += '_avg%i'%options.average
    if options.median > 0:
      elevation = scipy.ndimage.filters.median_filter(elevation,options.median)
      filters += '_median%i'%options.median

    np.savetxt(os.path.splitext(file)[0]+filters+os.path.splitext(file)[1],elevation,header='\n'.join(header))

