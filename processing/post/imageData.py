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
Generate PNG image from data in given column (or 2D data of overall table).

""", version = scriptID)

parser.add_option('-l','--label', dest='label', type='string',
                                  help='column containing data [all])')
parser.add_option('-r','--range', dest='range', type='float', nargs=2,
                                  help='data range (min max) [auto]')
parser.add_option('-d','--dimension', dest='dimension', type='int', nargs=2,
                                  help='data dimension (width height) [native]')
parser.add_option('--abs',        dest='abs', action='store_true',
                                  help='magnitude of values')
parser.add_option('--log',        dest='log', action='store_true',
                                  help='log10 of values')
parser.add_option('--fliplr',     dest='flipLR', action='store_true',
                                  help='flip around vertical axis')
parser.add_option('--flipud',     dest='flipUD', action='store_true',
                                  help='flip around horizontal axis')
parser.add_option('--color',      dest='color', type='string',
                                  help='color scheme')
parser.add_option('--invert',     dest='invert', action='store_true',
                                  help='invert color scheme')
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
                    range = [0.0,0.0],
                    dimension = [],
                    abs = False,
                    log = False,
                    flipLR = False,
                    flipUD = False,
                    color = "gray",
                    invert = False,
                    crop = [0,0,0,0],
                    pixelsize = 1,
                    pixelsizex = 1,
                    pixelsizey = 1,
                    show = False,
                   )

(options,filenames) = parser.parse_args()

if options.pixelsize > 1: (options.pixelsizex,options.pixelsizey) = [options.pixelsize]*2

# --- color palette ---------------------------------------------------------------------------------

theMap = damask.Colormap(predefined=options.color)
if options.invert: theMap = theMap.invert()
theColors = np.uint8(np.array(theMap.export(format='list',steps=256))*255)

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
            'output':open(os.path.splitext(name)[0]+\
                         ('_%s'%(options.label) if options.label != None else '')+\
                         '.png','w'),
            'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

  table = damask.ASCIItable(file['input'],file['output'],
                            buffered = False,                                                       # make unbuffered ASCII_table
                            labels = options.label != None)                                         # no labels when taking 2D dataset
  table.head_read()                                                                                 # read ASCII header info

# --------------- figure out column to process -----------------------------------------------------

  column = None
  if options.label != None:
    if options.label in table.labels:
      column = table.labels.index(options.label)                                                    # remember columns of requested data
    else:
      file['croak'].write('column %s not found...\n'%options.label)
      continue

# ------------------------------------------ process data ------------------------------------------

  table.data_readArray([column] if column != None else [])

  # convert data to values between 0 and 1 and arrange according to given options
  if options.dimension != []: table.data = table.data.reshape(options.dimension[1],options.dimension[0])
  if options.abs:             table.data = np.abs(table.data)
  if options.log:             table.data = np.log10(table.data)
  if options.flipLR:          table.data = np.fliplr(table.data)
  if options.flipUD:          table.data = np.flipud(table.data)
  if np.all(np.array(options.range) == 0.0): options.range = [table.data.min(),table.data.max()]

  table.data =         (table.data - min(options.range)) / \
               (max(options.range) - min(options.range))

  table.data = np.clip(table.data,0.0,1.0).\
                  repeat(options.pixelsizex,axis=1).\
                  repeat(options.pixelsizey,axis=0)

  (height,width) = table.data.shape

  im = Image.fromarray(theColors[np.array(255*table.data,dtype='i')], 'RGB').\
             crop((       options.crop[0],
                          options.crop[2],
                   width -options.crop[1],
                   height-options.crop[3]))

# ------------------------------------------ output result -----------------------------------------
  im.save(file['output'],format = "PNG")
  if options.show: im.show()

  table.input_close()                                                                               # close input ASCII table
  table.output_close()                                                                              # close output
