#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from optparse import OptionParser
from PIL import Image, ImageDraw
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Generate PNG image from data in given column vector containing RGB tuples.

""", version = scriptID)

parser.add_option('-l','--label', dest='label', type='string',
                                  help='column containing RGB triplet')
parser.add_option('-d','--dimension', dest='dimension', type='int', nargs=2,
                                  help='data dimension (width height)')
parser.add_option('--fliplr',     dest='flipLR', action='store_true',
                                  help='flip around vertical axis')
parser.add_option('--flipud',     dest='flipUD', action='store_true',
                                  help='flip around horizontal axis')
parser.add_option('--crop',       dest='crop', type='int', nargs=4, metavar='LEFT RIGHT TOP BOTTOM',
                                  help='pixels cropped on left, right, top, bottom')
parser.add_option('--show',       dest='show', action='store_true',
                                  help='show resulting image')
parser.add_option('-N','--pixelsize', dest='pixelsize', type='int',
                                  help='pixel per data point')
parser.add_option('-x','--pixelsizex', dest='pixelsizex', type='int',
                                  help='pixel per data point along x')
parser.add_option('-y','--pixelsizey', dest='pixelsizey', type='int',
                                  help='pixel per data point along y')

parser.set_defaults(label = None,
                    dimension = [],
                    flipLR = False,
                    flipUD = False,
                    crop = [0,0,0,0],
                    pixelsize = 1,
                    pixelsizex = 1,
                    pixelsizey = 1,
                    show = False,
                   )

(options,filenames) = parser.parse_args()

if options.dimension == []: parser.error('dimension of data array missing')
if options.pixelsize > 1: (options.pixelsizex,options.pixelsizey) = [options.pixelsize]*2

# --- loop over input files -------------------------------------------------------------------------
if filenames == []:
  filenames = ['STDIN']

for name in filenames:
  if name == 'STDIN':
    file = {'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m\n')
  else:
    if not os.path.exists(name): continue
    file = {'name':name,
            'input':open(name),
            'output':open(os.path.splitext(name)[0]+ \
                         ('' if options.label == None else '_'+options.label)+ \
                         '.png','w'),
            'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],
                            buffered = False,                                                       # make unbuffered ASCII_table
                            labels = options.label != None)                                         # no labels when taking 2D dataset
  table.head_read()                                                                                 # read ASCII header info

# ------------------------------------------ process data ------------------------------------------

  missing_labels = table.data_readArray(["%i_%s"%(i,options.label) for i in [1,2,3]])

  if len(missing_labels) > 0:
    file['croak'].write('column%s %s not found...\n'%('s' if len(missing_labels) > 1 else '',', '.join(missing_labels)))
    table.close(dismiss = True)                                                                     # close ASCII table file handles and delete output file
    continue

  # convert data to shape and arrange according to given options
  if options.dimension != []: table.data = table.data.reshape(options.dimension[1],options.dimension[0],3)
  if options.flipLR:          table.data = np.fliplr(table.data)
  if options.flipUD:          table.data = np.flipud(table.data)

  table.data = table.data.\
                          repeat(options.pixelsizex,axis=1).\
                          repeat(options.pixelsizey,axis=0)

  table.data *= 1. if np.any(table.data > 1.0) else 255.0                                          # ensure 8 bit data range

  (height,width,bands) = table.data.shape

  im = Image.fromarray(table.data.astype('uint8'), 'RGB').\
             crop((       options.crop[0],
                          options.crop[2],
                   width -options.crop[1],
                   height-options.crop[3]))

# ------------------------------------------ output result -----------------------------------------

  im.save(file['output'],format = "PNG")
  if options.show: im.show()

  table.close()                                                                                    # close ASCII table file handles
