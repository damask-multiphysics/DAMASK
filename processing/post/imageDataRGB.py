#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from optparse import OptionParser
from PIL import Image
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Generate PNG image from data in given column vector containing RGB tuples.

""", version = scriptID)

parser.add_option('-l','--label',
                 dest = 'label',
                 type = 'string', metavar = 'string',
                 help = 'column containing RGB triplet')
parser.add_option('-d','--dimension',
                  dest = 'dimension',
                  type = 'int', nargs = 2, metavar = 'int int',
                  help = 'data dimension (width height)')
parser.add_option('--fliplr',
                  dest = 'flipLR',
                  action = 'store_true',
                  help = 'flip around vertical axis')
parser.add_option('--flipud',
                  dest = 'flipUD',
                  action = 'store_true',
                  help = 'flip around horizontal axis')
parser.add_option('--crop',
                  dest = 'crop',
                  type = 'int', nargs = 4, metavar = ' '.join(['int']*4),
                  help = 'pixels cropped on left, right, top, bottom')
parser.add_option('-N','--pixelsize',
                  dest = 'pixelsize',
                  type = 'int', metavar = 'int',
                  help = 'pixels per data point')
parser.add_option('-x','--pixelsizex',
                  dest = 'pixelsizex',
                  type = 'int', metavar = 'int',
                  help = 'pixels per data point along x')
parser.add_option('-y','--pixelsizey',
                  dest = 'pixelsizey',
                  type = 'int', metavar = 'int',
                  help = 'pixels per data point along y')
parser.add_option('--show',
                  dest = 'show',
                  action = 'store_true',
                  help = 'show resulting image')

parser.set_defaults(label = None,
                    dimension = [],
                    flipLR = False,
                    flipUD = False,
                    crop = [0,0,0,0],
                    pixelsize  = 1,
                    pixelsizex = 1,
                    pixelsizey = 1,
                    show = False,
                   )

(options,filenames) = parser.parse_args()

if options.dimension == []: parser.error('dimension of data array missing')
if options.pixelsize > 1: (options.pixelsizex,options.pixelsizey) = [options.pixelsize]*2

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              buffered = False,
                              labeled = options.label != None,
                              readonly = True)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ process data ------------------------------------------

  errors = []
  
  missing_labels = table.data_readArray(options.label)
  if len(missing_labels) > 0:
    errors.append('column{} {} not found'.format('s' if len(missing_labels) > 1 else '',
                                                 ', '.join(missing_labels)))
  if table.label_dimension(options.label) != 3:
    errors.append('column {} has wrong dimension'.format(options.label))

  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss = True)                                                                     # close ASCII table file handles and delete output file
    continue

  # convert data to shape and arrange according to given options
  if options.dimension != []: table.data = table.data.reshape(options.dimension[1],options.dimension[0],3)
  if options.flipLR:          table.data = np.fliplr(table.data)
  if options.flipUD:          table.data = np.flipud(table.data)

  table.data = table.data.repeat(options.pixelsizex,axis=1).\
                          repeat(options.pixelsizey,axis=0)

  table.data *= 1. if np.any(table.data > 1.0) else 255.0                                          # ensure 8 bit data range

  (height,width,bands) = table.data.shape
  damask.util.croak('image dimension: {0} x {1}'.format(width,height))

  im = Image.fromarray(table.data.astype('uint8'), 'RGB').\
             crop((       options.crop[0],
                          options.crop[2],
                   width -options.crop[1],
                   height-options.crop[3]))

# ------------------------------------------ output result -----------------------------------------

  im.save(os.path.splitext(name)[0]+ \
          ('_'+options.label if options.label else '')+ \
          '.png' if name else sys.stdout,
          format = "PNG")

  table.close()                                                                                     # close ASCII table
  if options.show: im.show()
