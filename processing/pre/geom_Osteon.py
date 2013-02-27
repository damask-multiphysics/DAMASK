#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math,numpy,time
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP


# -----------------------------
class extendedOption(Option):
# -----------------------------
# used for definition of new option parser action 'extend', which enables to take multiple option arguments
# taken from online tutorial http://docs.python.org/library/optparse.html
    
    ACTIONS = Option.ACTIONS + ("extend",)
    STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
    TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
    ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

    def take_action(self, action, dest, opt, value, values, parser):
        if action == "extend":
            lvalue = value.split(",")
            values.ensure_value(dest, []).extend(lvalue)
        else:
            Option.take_action(self, action, dest, opt, value, values, parser)


# ----------------------- MAIN -------------------------------



parser = OptionParser(option_class=extendedOption, usage='%prog', description = """
Generate a geometry file of an osteon enclosing the Harvesian canal and separated by interstitial tissue.
The osteon phase is lamellar with a twisted plywood structure.
Its fiber orientation is oscillating by +/- amplitude within one period.
""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-r', '--resolution', dest='resolution', type='int', nargs=2, \
                                    help='resolution (a,b) of grid')
parser.add_option('-d', '--dimension', dest='dimension', type='float', nargs=2, \
                                    help='physical dimension (x,y) of periodic patch')
parser.add_option('-c', '--canal',  dest='canal', type='float', \
                                    help='Haversian canal radius')
parser.add_option('-o', '--osteon', dest='osteon', type='float', \
                                    help='osteon radius (horizontal)')
parser.add_option('-l', '--lamella', dest='period', type='float', \
                                    help='lamella width')
parser.add_option('-a', '--amplitude', dest='amplitude', type='float', \
                                    help='amplitude of twisted plywood wiggle in deg')
parser.add_option(      '--aspect', dest='aspect', type='float', \
                                    help='osteon aspect ratio (vert/horiz)')
parser.add_option('-w', '--omega',  dest='omega', type='float', \
                                    help='rotation angle (around normal) of osteon')
parser.add_option('--homogenization', dest='homogenization', type='int', \
                                    help='homogenization index to be used')
parser.add_option('--crystallite',  dest='crystallite', type='int', \
                                    help='crystallite index to be used')
parser.add_option('--configuration', dest='config', action='store_true', \
                                     help='output material configuration')
parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                                    help='use two-dimensional geom data arrangement')

parser.set_defaults(canal = 25e-6)
parser.set_defaults(osteon = 100e-6)
parser.set_defaults(aspect = 1.0)
parser.set_defaults(omega = 0.0)
parser.set_defaults(period = 5e-6)
parser.set_defaults(amplitude = 60)
parser.set_defaults(dimension = numpy.array([300e-6,300e-6],'d'))
parser.set_defaults(resolution = numpy.array([512,512],'i'))
parser.set_defaults(homogenization = 1)
parser.set_defaults(crystallite = 1)
parser.set_defaults(config = False)
parser.set_defaults(twoD = False)

(options, args) = parser.parse_args()

# ------------------------------------------ setup file handles ---------------------------------------  

file = {'name':'STDIN',
        'input':sys.stdin,
        'output':sys.stdout,
        'croak':sys.stderr,
       }

if numpy.any(options.resolution < 2):
  file['croak'].write('resolution too low...\n')
  sys.exit()

options.omega  *= math.pi/180.0                                           # rescale ro radians
rotation = numpy.array([[ math.cos(options.omega),math.sin(options.omega),],
                        [-math.sin(options.omega),math.cos(options.omega),]],'d')

box = numpy.dot(numpy.array([[options.canal,0.],[0.,options.aspect*options.canal]]).transpose(),rotation)
sys.stderr.write("bounding box: %s\n"%(numpy.sqrt(numpy.sum(box*box,0))))

info = {'grains': 0,
        'resolution': numpy.ones(3,'i'),
        'dimension':  numpy.ones(3,'d'),
        'origin':     numpy.zeros(3,'d'),
        'homogenization': options.homogenization,
       }


info['resolution'][:2] = options.resolution
info['dimension'][:2]  = options.dimension
info['dimension'][2]   = min(info['dimension'][0]/info['resolution'][0],info['dimension'][1]/info['resolution'][1])
info['origin']         = -info['dimension']/2.0

X0 = info['dimension'][0]/info['resolution'][0]*\
     (numpy.tile(numpy.arange(info['resolution'][0]),(info['resolution'][1],1))             - info['resolution'][0]/2 + 0.5)
Y0 = info['dimension'][1]/info['resolution'][1]*\
     (numpy.tile(numpy.arange(info['resolution'][1]),(info['resolution'][0],1)).transpose() - info['resolution'][1]/2 + 0.5)

X = X0*rotation[0,0] + Y0*rotation[0,1]                                   # rotate by omega
Y = X0*rotation[1,0] + Y0*rotation[1,1]                                   # rotate by omega

radius = numpy.sqrt(X*X + Y*Y/options.aspect/options.aspect)
alpha = numpy.degrees(numpy.arctan2(Y/options.aspect,X))
beta = options.amplitude*numpy.sin(2.0*math.pi*(radius-options.canal)/options.period)

microstructure = numpy.where(radius < float(options.canal),1,0) + numpy.where(radius > float(options.osteon),2,0)

info['grains'] = 3
alphaOfGrain = numpy.zeros(info['resolution'][0]*info['resolution'][1],'d')
betaOfGrain  = numpy.zeros(info['resolution'][0]*info['resolution'][1],'d')
for y in xrange(info['resolution'][1]):
  for x in xrange(info['resolution'][0]):
    if microstructure[y,x] == 0:
      microstructure[y,x] = info['grains']
      alphaOfGrain[info['grains']] = alpha[y,x]
      betaOfGrain[ info['grains']] = beta[y,x]
      info['grains'] += 1


# -------------------------------------- switch according to task ----------------------------------

formatwidth = 1+int(math.floor(math.log10(info['grains']-1)))

if options.config:
  file['output'].write('<microstructure>\n')
  file['output'].write('\n[canal]\n' + \
                       'crystallite %i\n'%options.crystallite + \
                       '(constituent)\tphase 1\ttexture 1\tfraction 1.0\n')
  file['output'].write('\n[interstitial]\n' + \
                       'crystallite %i\n'%options.crystallite + \
                       '(constituent)\tphase 2\ttexture 2\tfraction 1.0\n')
  for i in xrange(3,info['grains']):
    file['output'].write('\n[Grain%s]\n'%(str(i).zfill(formatwidth)) + \
                         'crystallite %i\n'%options.crystallite + \
                         '(constituent)\tphase 3\ttexture %s\tfraction 1.0\n'%(str(i).rjust(formatwidth)))

  file['output'].write('\n<texture>\n')
  file['output'].write('\n[canal]\n')
  file['output'].write('\n[interstitial]\n')
  for i in xrange(3,info['grains']):
    file['output'].write('\n[Grain%s]\n'%(str(i).zfill(formatwidth)) + \
                         '(gauss)\tphi1 %g\tPhi %g\tphi2 0\tscatter 0.0\tfraction 1.0\n'%(\
                                  alphaOfGrain[i],\
                                  betaOfGrain[i]))
  
else:

  file['output'].write("4 header\n" + \
                       "resolution\ta %i\tb %i\tc %i\n"%(info['resolution'][0],info['resolution'][1],info['resolution'][2]) + \
                       "dimension\tx %g\ty %g\tz %g\n"%(info['dimension'][0],info['dimension'][1],info['dimension'][2]) + \
                       "origin\tx %g\ty %g\tz %g\n"%(info['origin'][0],info['origin'][1],info['origin'][2]) + \
                       "homogenization 1\n"
                       )
  
  for y in xrange(info['resolution'][1]):
    for x in xrange(info['resolution'][0]):
      file['output'].write(\
          str(microstructure[y,x]).rjust(formatwidth) + \
          {True:' ',False:'\n'}[options.twoD] )
    file['output'].write({True:'\n',False:''}[options.twoD])
  

# ------------------------------------------ output finalization ---------------------------------------  

file['output'].close()
