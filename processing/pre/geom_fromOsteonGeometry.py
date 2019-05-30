#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser

import numpy as np

import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [geomfile]', description = """
Generate description of an osteon enclosing the Harvesian canal and separated by interstitial tissue.
The osteon phase is lamellar with a twisted plywood structure.
Its fiber orientation is oscillating by +/- amplitude within one period.

""", version = scriptID)


parser.add_option('-g', '--grid',
                  dest='grid', type='int',
                  nargs=2, metavar = 'int int',
                  help='a,b grid of hexahedral box [%default]')
parser.add_option('-s', '--size',
                  dest='size',
                  type='float', nargs=2, metavar = 'float float',
                  help='x,y size of hexahedral box [%default]')
parser.add_option('-c', '--canal',
                  dest='canal',
                  type='float', metavar = 'float',
                  help='Haversian canal radius [%default]')
parser.add_option('-o', '--osteon',
                  dest='osteon',
                  type='float', metavar = 'float',
                  help='horizontal osteon radius [%default]')
parser.add_option('-l', '--lamella',
                  dest='period',
                  type='float', metavar = 'float',
                  help='lamella width [%default]')
parser.add_option('-a', '--amplitude',
                  dest='amplitude',
                  type='float', metavar = 'float',
                  help='amplitude of twisted plywood wiggle in deg [%default]')
parser.add_option(      '--aspect',
                  dest='aspect',
                  type='float', metavar = 'float',
                  help='vertical/horizontal osteon aspect ratio [%default]')
parser.add_option('-w', '--omega',
                  dest='omega',
                  type='float', metavar = 'float',
                  help='rotation angle around normal of osteon [%default]')
parser.add_option(      '--homogenization',
                  dest='homogenization',
                  type='int', metavar = 'int',
                  help='homogenization index to be used [%default]')

parser.set_defaults(canal = 25e-6,
                    osteon = 100e-6,
                    aspect = 1.0,
                    omega = 0.0,
                    period = 5e-6,
                    amplitude = 60,
                    size = (300e-6,300e-6),
                    grid = (512,512),
                    homogenization = 1)

(options,filename) = parser.parse_args()


name = None if filename == [] else filename[0]
damask.util.report(scriptName,name)

omega = np.deg2rad(options.omega)
rotation = np.array([[ np.cos(omega),np.sin(omega),],
                     [-np.sin(omega),np.cos(omega),]])

grid = np.array(options.grid,'i')
size = np.array(options.size,'d')

X0,Y0 = np.meshgrid(size[0]/grid[0] * (np.arange(grid[0]) - grid[0]/2 + 0.5),
                    size[1]/grid[0] * (np.arange(grid[1]) - grid[1]/2 + 0.5), indexing='ij')

X = X0*rotation[0,0] + Y0*rotation[0,1]                                                             # rotate by omega
Y = X0*rotation[1,0] + Y0*rotation[1,1]                                                             # rotate by omega

radius = np.sqrt(X*X + Y*Y/options.aspect**2.0)
alpha = np.degrees(np.arctan2(Y/options.aspect,X))
beta  = options.amplitude*np.sin(2.0*np.pi*(radius-options.canal)/options.period)

microstructure = np.where(radius < float(options.canal), 1,0) \
               + np.where(radius > float(options.osteon),2,0)

# extend to 3D
size = np.append(size,np.min(size/grid))
grid = np.append(grid,1)
microstructure = microstructure.reshape(microstructure.shape+(1,))

Alpha = np.zeros(grid[0]*grid[1],'d')
Beta  = np.zeros(grid[0]*grid[1],'d')

i = 3
for y in range(grid[1]):
  for x in range(grid[0]):
    if microstructure[x,y] == 0:
      microstructure[x,y] = i
      Alpha[i] = alpha[x,y]
      Beta [i] = beta [x,y]
      i+=1

config_header = ['<texture>',
                 '[canal]',
                 '[interstitial]'
                ]
for i in range(3,np.max(microstructure)):
  config_header += ['[Point{}]'.format(i-2),
                    '(gauss)\tphi1 {:.2f}\tPhi {:.2f}\tphi2 0'.format(Alpha[i],Beta[i])
                   ]

config_header = ['<microstructure>',
                 '[canal]',
                 'crystallite 1',
                 '(constituent)\tphase 1\ttexture 1\tfraction 1.0',
                 '[interstitial]',
                 'crystallite 1',
                 '(constituent)\tphase 2\ttexture 2\tfraction 1.0'
                ]
for i in range(3,np.max(microstructure)):
  config_header += ['[Point{}]'.format(i-2),
                    'crystallite 1',
                    '(constituent)\tphase 3\ttexture {}\tfraction 1.0'.format(i)
                   ]

header = [scriptID + ' ' + ' '.join(sys.argv[1:])]\
       + config_header
geom = damask.Geom(microstructure.reshape(grid),
                   size,-size/2,
                   homogenization=options.homogenization,comments=header)
damask.util.croak(geom)

if name is None:
  sys.stdout.write(str(geom.show()))
else:
  geom.to_file(name)
