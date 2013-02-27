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

parser.add_option('-c', '--canal',  dest='canal', type='float', \
                                    help='Haversian canal radius')
parser.add_option('-o', '--osteon', dest='osteon', type='float', \
                                    help='Osteon radius (horizontal)')
parser.add_option('-r', '--aspect', dest='aspect', type='float', \
                                    help='Osteon aspect ratio (vert/horiz)')
parser.add_option('-w', '--omega',  dest='omega', type='float', \
                                    help='rotation angle (around normal) of osteon')
parser.add_option('-s', '--size',   dest='size', type='float', \
                                    help='box size (horizontal)')
parser.add_option('-m', '--margin', dest='margin', type='float', \
                                    help='margin width')
parser.add_option('-N',  '--resolution', dest='resolution', type='int', \
                                    help='box size in pixels (horizontal)')
parser.add_option('-a', '--amplitude', dest='amplitude', type='float', \
                                    help='Amplitude of twisted plywood wiggle in deg')
parser.add_option('-p', '--period', dest='period', type='float', \
                                    help='lamella width')
parser.add_option('--homogenization', dest='homogenization', type='int', \
                                    help='homogenization index to be used')
parser.add_option('--crystallite',  dest='crystallite', type='int', \
                                    help='crystallite index to be used')
parser.add_option('--configuration', dest='config', action='store_true', \
                                     help='output material configuration')
parser.add_option('-2', '--twodimensional', dest='twoD', action='store_true', \
                                    help='output geom file with two-dimensional data arrangement')

parser.set_defaults(canal = 25e-6)
parser.set_defaults(osteon = 75e-6)
parser.set_defaults(aspect = 1.0)
parser.set_defaults(omega = 0.0)
parser.set_defaults(period = 5e-6)
parser.set_defaults(amplitude = 60)
parser.set_defaults(size = 300e-6)
parser.set_defaults(margin = 0.0)
parser.set_defaults(resolution = 256)
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

if (options.resolution < 2):
  file['croak'].write('resolution too low...\n')
  sys.exit()

if (options.size < options.canal):
  file['croak'].write('canal larger than box...\n')
  sys.exit()

if (options.osteon < options.canal):
  file['croak'].write('canal larger than osteon...\n')
  sys.exit()

info = {'grains': 0,
        'resolution': numpy.zeros(3,'i'),
        'dimension':  numpy.ones(3,'d'),
        'origin':     numpy.ones(3,'d'),
        'homogenization': options.homogenization,
       }

if options.margin > 0.0:
  info['resolution'][0] = options.resolution
  info['resolution'][1] = round(options.resolution*(options.aspect * options.osteon + options.margin) / \
                                                                    (options.osteon + options.margin) )
  info['dimension'][0]  =                  options.osteon + options.margin
  info['dimension'][1]  = options.aspect * options.osteon + options.margin
else:
  info['resolution'][0] = options.resolution
  info['resolution'][1] = options.resolution
  info['dimension'][0]  = options.size
  info['dimension'][1]  = options.size

options.canal  *= info['resolution'][0]/info['dimension'][0]              # rescale to pixel dimension
options.osteon *= info['resolution'][0]/info['dimension'][0]              # rescale to pixel dimension
options.period *= info['resolution'][0]/info['dimension'][0]              # rescale to pixel dimension
options.omega  *= math.pi/180.0                                           # rescale ro radians

rotation = numpy.array([[ math.cos(options.omega),math.sin(options.omega),],
                        [-math.sin(options.omega),math.cos(options.omega),]],'d')

X0 = numpy.tile(range(info['resolution'][0]),(info['resolution'][0],1))             - info['resolution'][0]/2 + 0.5
Y0 = numpy.tile(range(info['resolution'][1]),(info['resolution'][1],1)).transpose() - info['resolution'][1]/2 + 0.5

X = X0*rotation[0,0] + Y0*rotation[1,0]                                   # rotate by omega
Y = X0*rotation[0,1] + Y0*rotation[1,1]                                   # rotate by omega

radius = numpy.sqrt(X*X/options.aspect/options.aspect + Y*Y)
alpha = numpy.degrees(numpy.arctan2(X/options.aspect,Y))
beta = options.amplitude*numpy.sin(math.pi*(radius-options.canal)/options.period)


microstructure = numpy.where(radius < float(options.canal),1,0) + numpy.where(radius > float(options.osteon),2,0)
sys.stderr.write('micro %f\n'%(time.clock()))

info['grains'] = 3
alphaOfGrain = numpy.zeros(info['resolution'][0]*info['resolution'][1],'d')
betaOfGrain  = numpy.zeros(info['resolution'][0]*info['resolution'][1],'d')
for y in xrange(info['resolution'][1]):
  for x in xrange(info['resolution'][0]):
    if microstructure[x,y] == 0:
      microstructure[x,y] = info['grains']
      alphaOfGrain[info['grains']] = alpha[x,y]
      betaOfGrain[ info['grains']] = beta[x,y]
      info['grains'] += 1


sys.stderr.write('micro assemble %f\n'%(time.clock()))


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
                       "resolution\ta %i\tb %i\tc %i\n"%(info['resolution'][0],info['resolution'][1],1,) + \
                       "dimension\tx %g\ty %g\tz %g\n"%(info['dimension'][0],info['dimension'][1],info['dimension'][0]/info['resolution'][0],) + \
                       "origin\tx 0\ty 0\tz 0\n" + \
                       "homogenization 1\n"
                       )
  
  for y in xrange(info['resolution'][1]):
    for x in xrange(info['resolution'][0]):
      file['output'].write(\
          str(microstructure[x,y]).rjust(formatwidth) + \
          {True:' ',False:'\n'}[options.twoD] )
    file['output'].write({True:'\n',False:''}[options.twoD])
  

# ------------------------------------------ output finalization ---------------------------------------  

file['output'].close()
