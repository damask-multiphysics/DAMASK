#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string
import numpy as np
from optparse import OptionParser
from PIL import Image, ImageDraw
import damask

scriptID   = string.replace('$Id: imageDataDeformed.py 4188 2015-05-20 23:21:35Z p.eisenlohr $','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Generate PNG image from scalar data on grid deformed by (periodic) deformation gradient.

""", version = scriptID)

parser.add_option('-l','--label', dest='label', type='string', metavar='string',
                                  help='column containing data)')
parser.add_option('-r','--range', dest='range', type='float', nargs=2, metavar='float float',
                                  help='data range (min max) [auto]')
parser.add_option('--color',      dest='color', type='string', metavar='string',
                                  help='color scheme')
parser.add_option('--invert',     dest='invert', action='store_true',
                                  help='invert color scheme')
parser.add_option('--abs',        dest='abs', action='store_true',
                                  help='magnitude of values')
parser.add_option('--log',        dest='log', action='store_true',
                                  help='log10 of values')
parser.add_option('-d','--dimension', dest='dimension', type='int', nargs=3, metavar=' '.join(['int']*3),
                                  help='data dimension (x/y/z)')
parser.add_option('-s','--size',  dest='size', type='float', nargs=3, metavar=' '.join(['float']*3),
                                  help='box size (x/y/z)')
parser.add_option('-f','--defgrad',   dest='defgrad', metavar='string',
                                  help='column label of deformation gradient [%default]')
parser.add_option('--scaling',    dest='scaling', type='float', nargs=3, metavar = ' '.join(['float']*3),
                                  help='x/y/z scaling of displacment fluctuation [%default]')
parser.add_option('-z','--layer', dest='z', type='int', metavar='int',
                                  help='index of z plane to plot [%default]')
parser.add_option('--fliplr',     dest='flipLR', action='store_true',
                                  help='flip around vertical axis')
parser.add_option('--flipud',     dest='flipUD', action='store_true',
                                  help='flip around horizontal axis')
parser.add_option('--crop',       dest='crop', type='int', nargs=4, metavar=' '.join(['int']*3),
                                  help='pixels cropped on left, right, top, bottom')
parser.add_option('--show',       dest='show', action='store_true',
                                  help='show resulting image')
parser.add_option('-N','--pixelsize', dest='pixelsize', type='int', metavar='int',
                                  help='pixels per cell edge')

parser.set_defaults(label = None,
                    range = [0.0,0.0],
                    dimension = [],
                    size = [],
                    z = 1,
                    abs = False,
                    log = False,
                    defgrad = 'f',
                    scaling = [1.,1.,1.],
                    flipLR = False,
                    flipUD = False,
                    color = "gray",
                    invert = False,
                    crop = [0,0,0,0],
                    pixelsize = 1,
                    show = False,
                   )

(options,filenames) = parser.parse_args()

options.size      = np.array(options.size)
options.dimension = np.array(options.dimension)
options.range     = np.array(options.range)

if options.z > 0: options.z -= 1                                                                    # adjust to 0-based indexing

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
            'output':open(os.path.splitext(name)[0]+ \
                         ('' if options.label == None else '_'+options.label)+ \
                         '.png','w'),
            'croak':sys.stderr}
    file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')
  
  table = damask.ASCIItable(file['input'],file['output'],
                            buffered = False,                                                       # make unbuffered ASCII_table
                            labels = options.label != None)                                         # no labels when taking 2D dataset
  table.head_read()                                                                                 # read ASCII header info
  
# --------------- figure out columns to process  ---------------------------------------------------
  
  errors = []
  if table.label_dimension(options.label) != 1:
    errors.append('no scalar data (%s) found...'%options.label)
  if table.label_dimension(options.defgrad) != 9:
    errors.append('no deformation gradient tensor (1..9_%s) found...'%options.defgrad)
  
  if errors != []:
    file['croak'].write('\n'.join(errors)+'\n')
    table.close(dismiss = True)
    continue
  
  table.data_readArray([options.label,options.defgrad])
  
  F =    table.data[:,1:10].transpose().reshape([3,3]+list(options.dimension),order='F')
  data = table.data[:,0   ].transpose().reshape(      list(options.dimension),order='F')

  if options.abs:  data = np.abs(data)
  if options.log:  data = np.log10(data)
  if np.all(options.range == 0.0): options.range = np.array([data.min(),data.max()])
  elif options.log:                options.range = np.log10(options.range)
  
  data = (       data         - options.range.min()) / \
         (options.range.max() - options.range.min())                                                  # data scaled to fraction of range
  
  data = np.clip(data,0.0,1.0)                                                                        # cut off outliers (should be none)
  
# ---------------- calculate coordinates -----------------------------------------------------------
  
  Favg      = damask.core.math.tensorAvg(F)
  centroids = damask.core.mesh.deformedCoordsFFT(options.size,F,Favg,options.scaling)
  nodes     = damask.core.mesh.nodesAroundCentres(options.size,Favg,centroids)
  
  boundingBox = np.array([ \
                 [np.amin(nodes[0,:,:,options.z]),np.amin(nodes[1,:,:,options.z]),np.amin(nodes[2,:,:,options.z])],
                 [np.amax(nodes[0,:,:,options.z]),np.amax(nodes[1,:,:,options.z]),np.amax(nodes[2,:,:,options.z])],
                ])                                                                                  # find x-y bounding box for given z layer
  
  nodes -= boundingBox[0].repeat(np.prod(options.dimension+1)).reshape([3]+list(options.dimension+1))
  nodes *= (options.pixelsize*options.dimension/options.size).repeat(np.prod(options.dimension+1)).reshape([3]+list(options.dimension+1))
  imagesize = (options.pixelsize*(boundingBox[1]-boundingBox[0])*options.dimension\
                                                                /options.size)[:2].astype('i')      # determine image size from number of cells in overall bounding box
  im = Image.new('RGBA',imagesize)
  draw = ImageDraw.Draw(im)
  
  for y in xrange(options.dimension[1]):
    for x in xrange(options.dimension[0]):
      draw.polygon([nodes[0,x  ,y  ,options.z],
                    nodes[1,x  ,y  ,options.z],
                    nodes[0,x+1,y  ,options.z],
                    nodes[1,x+1,y  ,options.z],
                    nodes[0,x+1,y+1,options.z],
                    nodes[1,x+1,y+1,options.z],
                    nodes[0,x  ,y+1,options.z],
                    nodes[1,x  ,y+1,options.z],
                   ],
                   fill =    tuple(theColors[int(255*data[x,y,options.z])]),
                   outline = None)
  
#  if options.flipLR:          table.data = np.fliplr(table.data)
#  if options.flipUD:          table.data = np.flipud(table.data)


#  (height,width,bands) = table.data.shape

#  im = Image.fromarray(table.data.astype('uint8'), 'RGB').\
#             crop((       options.crop[0],
#                          options.crop[2],
#                   width -options.crop[1],
#                   height-options.crop[3]))

# ------------------------------------------ output result -----------------------------------------
  
  im.save(file['output'],format = "PNG")
  if options.show: im.show()
  
  table.close()                                                                                    # close ASCII table file handles
