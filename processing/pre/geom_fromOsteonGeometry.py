#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [options]', description = """
Generate a geometry file of an osteon enclosing the Harvesian canal and separated by interstitial tissue.
The osteon phase is lamellar with a twisted plywood structure.
Its fiber orientation is oscillating by +/- amplitude within one period.

""", version = scriptID)


parser.add_option('-g', '--grid', dest='grid', type='int', nargs=2, metavar = 'int int',
                  help='a,b grid of hexahedral box [%default]')
parser.add_option('-s', '--size', dest='size', type='float', nargs=2, metavar = 'float float',
                  help='x,y size of hexahedral box [%default]')
parser.add_option('-c', '--canal',  dest='canal', type='float', metavar = 'float',
                  help='Haversian canal radius [%default]')
parser.add_option('-o', '--osteon', dest='osteon', type='float', metavar = 'float',
                  help='horizontal osteon radius [%default]')
parser.add_option('-l', '--lamella', dest='period', type='float', metavar = 'float',
                  help='lamella width [%default]')
parser.add_option('-a', '--amplitude', dest='amplitude', type='float', metavar = 'float',
                  help='amplitude of twisted plywood wiggle in deg [%default]')
parser.add_option(      '--aspect', dest='aspect', type='float', metavar = 'float',
                  help='vertical/horizontal osteon aspect ratio [%default]')
parser.add_option('-w', '--omega',  dest='omega', type='float', metavar = 'float',
                  help='rotation angle around normal of osteon [%default]')
parser.add_option('--homogenization', dest='homogenization', type='int', metavar = 'int',
                  help='homogenization index to be used [%default]')
parser.add_option('--crystallite',  dest='crystallite', type='int', metavar = 'int',
                  help='crystallite index to be used [%default]')

parser.set_defaults(canal = 25e-6,
                    osteon = 100e-6,
                    aspect = 1.0,
                    omega = 0.0,
                    period = 5e-6,
                    amplitude = 60,
                    size = (300e-6,300e-6),
                    grid = (512,512),
                    homogenization = 1,
                    crystallite = 1)

(options,filename) = parser.parse_args()

if np.any(options.grid < 2):
  parser('invalid grid a b c.')
if np.any(options.size <= 0.0):
  parser('invalid size x y z.')

# --- open input files ----------------------------------------------------------------------------

if filename == []: filename = [None]

table = damask.ASCIItable(outname = filename[0],
                          buffered = False)

damask.util.report(scriptName,filename[0])

options.omega  *= math.pi/180.0                                                                     # rescale ro radians
rotation = np.array([[ math.cos(options.omega),math.sin(options.omega),],
                     [-math.sin(options.omega),math.cos(options.omega),]],'d')

box = np.dot(np.array([[options.canal,0.],[0.,options.aspect*options.canal]]).transpose(),rotation)


info = {
        'grid':   np.ones(3,'i'),
        'size':   np.ones(3,'d'),
        'origin': np.zeros(3,'d'),
        'microstructures': 3,
        'homogenization':  options.homogenization,
       }

info['grid'][:2] = np.array(options.grid,'i')
info['size'][:2] = np.array(options.size,'d')
info['size'][2]  = min(info['size'][0]/info['grid'][0],info['size'][1]/info['grid'][1])
info['origin']   = -info['size']/2.0

X0 = info['size'][0]/info['grid'][0]*\
     (np.tile(np.arange(info['grid'][0]),(info['grid'][1],1))             - info['grid'][0]/2 + 0.5)
Y0 = info['size'][1]/info['grid'][1]*\
     (np.tile(np.arange(info['grid'][1]),(info['grid'][0],1)).transpose() - info['grid'][1]/2 + 0.5)

X = X0*rotation[0,0] + Y0*rotation[0,1]                                                             # rotate by omega
Y = X0*rotation[1,0] + Y0*rotation[1,1]                                                             # rotate by omega

radius = np.sqrt(X*X + Y*Y/options.aspect/options.aspect)
alpha = np.degrees(np.arctan2(Y/options.aspect,X))
beta = options.amplitude*np.sin(2.0*math.pi*(radius-options.canal)/options.period)

microstructure = np.where(radius < float(options.canal),1,0) + np.where(radius > float(options.osteon),2,0)

alphaOfGrain = np.zeros(info['grid'][0]*info['grid'][1],'d')
betaOfGrain  = np.zeros(info['grid'][0]*info['grid'][1],'d')
for y in xrange(info['grid'][1]):
  for x in xrange(info['grid'][0]):
    if microstructure[y,x] == 0:
      microstructure[y,x] = info['microstructures']
      alphaOfGrain[info['microstructures']] = alpha[y,x]
      betaOfGrain[ info['microstructures']] = beta[y,x]
      info['microstructures'] += 1

#--- report ---------------------------------------------------------------------------------------
damask.util.croak(['grid     a b c:  %s'%(' x '.join(map(str,info['grid']))),
                   'size     x y z:  %s'%(' x '.join(map(str,info['size']))),
                   'origin   x y z:  %s'%(' : '.join(map(str,info['origin']))),
                   'homogenization:  %i'%info['homogenization'],
                   'microstructures: %i'%info['microstructures']])
# -------------------------------------- switch according to task ----------------------------------
formatwidth = 1+int(math.floor(math.log10(info['microstructures']-1)))
header = [scriptID + ' ' + ' '.join(sys.argv[1:])]
header.append('<microstructure>')
header.append('[canal]')
header.append('crystallite %i'%options.crystallite)
header.append('(constituent)\tphase 1\ttexture 1\tfraction 1.0')
header.append('[interstitial]')
header.append('crystallite %i'%options.crystallite)
header.append('(constituent)\tphase 2\ttexture 2\tfraction 1.0')
for i in xrange(3,info['microstructures']):
  header.append('[Grain%s]'%(str(i).zfill(formatwidth)))
  header.append('crystallite %i'%options.crystallite)
  header.append('(constituent)\tphase 3\ttexture %s\tfraction 1.0'%(str(i).rjust(formatwidth)))

header.append('<texture>')
header.append('[canal]')
header.append('[interstitial]')
for i in xrange(3,info['microstructures']):
  header.append('[Grain%s]'%(str(i).zfill(formatwidth)))
  header.append('(gauss)\tphi1 %g\tPhi %g\tphi2 0\tscatter 0.0\tfraction 1.0'\
                                                          %(alphaOfGrain[i],betaOfGrain[i]))
header.append([
    "grid\ta {grid[0]}\tb {grid[1]}\tc {grid[2]}".format(grid=info['grid']),
    "size\tx {size[0]}\ty {size[1]}\tz {size[2]}".format(size=info['size']),
    "origin\tx {origin[0]}\ty {origin[1]}\tz {origin[2]}".format(origin=info['origin']),
    "homogenization\t{homog}".format(homog=info['homogenization']),
    "microstructures\t{microstructures}".format(microstructures=info['microstructures'])])
  
table.info_append(header)
table.head_write()
      
# --- write microstructure information ------------------------------------------------------------

table.data = microstructure.reshape(info['grid'][1]*info['grid'][2],info['grid'][0])
table.data_writeArray('%%%ii'%(formatwidth),delimiter=' ')
    
#--- output finalization --------------------------------------------------------------------------
table.close()  
