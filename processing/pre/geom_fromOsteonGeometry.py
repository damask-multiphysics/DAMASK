#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math
import numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog [options]', description = """
Generate a geometry file of an osteon enclosing the Harvesian canal and separated by interstitial tissue.
The osteon phase is lamellar with a twisted plywood structure.
Its fiber orientation is oscillating by +/- amplitude within one period.

""", version = scriptID)


parser.add_option('-g', '--grid', dest='grid', type='int', nargs=2, metavar = 'int int', \
                  help='a,b grid of hexahedral box %default')
parser.add_option('-s', '--size', dest='size', type='float', nargs=2, metavar = 'float float', \
                  help='x,y size of hexahedral box %default')
parser.add_option('-c', '--canal',  dest='canal', type='float', metavar = 'float', \
                  help='Haversian canal radius [%default]')
parser.add_option('-o', '--osteon', dest='osteon', type='float', metavar = 'float', \
                  help='osteon radius (horizontal) [%default]')
parser.add_option('-l', '--lamella', dest='period', type='float', metavar = 'float', \
                  help='lamella width [%default]')
parser.add_option('-a', '--amplitude', dest='amplitude', type='float', metavar = 'float', \
                  help='amplitude of twisted plywood wiggle in deg [%default]')
parser.add_option(      '--aspect', dest='aspect', type='float', metavar = 'float', \
                  help='osteon aspect ratio (vert/horiz) [%default]')
parser.add_option('-w', '--omega',  dest='omega', type='float', metavar = 'float', \
                  help='rotation angle (around normal) of osteon [%default]')
parser.add_option('--homogenization', dest='homogenization', type='int', metavar = 'int', \
                  help='homogenization index to be used [%default]')
parser.add_option('--crystallite',  dest='crystallite', type='int', metavar = 'int', \
                  help='crystallite index to be used [%default]')
parser.add_option('--configuration', dest='config', action='store_true', \
                  help='output material configuration [%default]')
parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                  help='use two-dimensional geom data arrangement [%default]')

parser.set_defaults(canal = 25e-6)
parser.set_defaults(osteon = 100e-6)
parser.set_defaults(aspect = 1.0)
parser.set_defaults(omega = 0.0)
parser.set_defaults(period = 5e-6)
parser.set_defaults(amplitude = 60)
parser.set_defaults(size = np.array([300e-6,300e-6],'d'))
parser.set_defaults(grid = np.array([512,512],'i'))
parser.set_defaults(homogenization = 1)
parser.set_defaults(crystallite = 1)
parser.set_defaults(config = False)
parser.set_defaults(twoD = False)

(options,filename) = parser.parse_args()

# ------------------------------------------ setup file handle -------------------------------------
if filename == []:
  file = {'output':sys.stdout, 'croak':sys.stderr}
else:
  file = {'output':open(filename[0],'w'), 'croak':sys.stderr}

if np.any(options.grid < 2):
  file['croak'].write('grid too small...\n')
  sys.exit()

if np.any(options.size <= 0.0):
  file['croak'].write('size too small...\n')
  sys.exit()

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

info['grid'][:2] = options.grid
info['size'][:2] = options.size
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
else: file['croak'].write('\033[1m'+scriptName+'\033[0m\n')
file['croak'].write('grid     a b c:  %s\n'%(' x '.join(map(str,info['grid']))) + \
                    'size     x y z:  %s\n'%(' x '.join(map(str,info['size']))) + \
                    'origin   x y z:  %s\n'%(' : '.join(map(str,info['origin']))) + \
                    'microstructures: %i\n'%info['microstructures'] + \
                    'homogenization:  %i\n'%info['homogenization'])
file['croak'].write("bounding box:    %s\n"%(np.sqrt(np.sum(box*box,0))))

if np.any(info['grid'] < 1):
  file['croak'].write('invalid grid a b c.\n')
  sys.exit()
if np.any(info['size'] <= 0.0):
  file['croak'].write('invalid size x y z.\n')
  sys.exit()

# -------------------------------------- switch according to task ----------------------------------
formatwidth = 1+int(math.floor(math.log10(info['microstructures']-1)))
if options.config:
  file['output'].write('<microstructure>\n')
  file['output'].write('\n[canal]\n' + \
                       'crystallite %i\n'%options.crystallite + \
                       '(constituent)\tphase 1\ttexture 1\tfraction 1.0\n')
  file['output'].write('\n[interstitial]\n' + \
                       'crystallite %i\n'%options.crystallite + \
                       '(constituent)\tphase 2\ttexture 2\tfraction 1.0\n')
  for i in xrange(3,info['microstructures']):
    file['output'].write('\n[Grain%s]\n'%(str(i).zfill(formatwidth)) + \
                         'crystallite %i\n'%options.crystallite + \
                         '(constituent)\tphase 3\ttexture %s\tfraction 1.0\n'%(str(i).rjust(formatwidth)))

  file['output'].write('\n<texture>\n')
  file['output'].write('\n[canal]\n')
  file['output'].write('\n[interstitial]\n')
  for i in xrange(3,info['microstructures']):
    file['output'].write('\n[Grain%s]\n'%(str(i).zfill(formatwidth)) + \
                         '(gauss)\tphi1 %g\tPhi %g\tphi2 0\tscatter 0.0\tfraction 1.0\n'%(\
                                  alphaOfGrain[i],\
                                  betaOfGrain[i]))
  
else:
  header = [scriptID + ' ' + ' '.join(sys.argv[1:])+'\n']
  header.append("grid\ta %i\tb %i\tc %i\n"%(info['grid'][0],info['grid'][1],info['grid'][2],))
  header.append("size\tx %f\ty %f\tz %f\n"%(info['size'][0],info['size'][1],info['size'][2],))
  header.append("origin\tx %f\ty %f\tz %f\n"%(info['origin'][0],info['origin'][1],info['origin'][2],))
  header.append("microstructures\t%i\n"%info['microstructures'])
  header.append("homogenization\t%i\n"%info['homogenization'])
  file['output'].write('%i\theader\n'%(len(header))+''.join(header))
  
  for y in xrange(info['grid'][1]):
    for x in xrange(info['grid'][0]):
      file['output'].write(\
          str(microstructure[y,x]).rjust(formatwidth) + \
          {True:' ',False:'\n'}[options.twoD] )
    file['output'].write({True:'\n',False:''}[options.twoD])
 
