#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math,numpy,random
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP

#------------------------------------------------------------------------------------------------
class extendedOption(Option):
#------------------------------------------------------------------------------------------------
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


#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
identifiers = {
        'grid': ['a','b','c'],
          }
mappings = {
        'grid': lambda x: int(x),
          }


parser = OptionParser(option_class=extendedOption, usage='%prog [options]', description = """
Distribute given number of points randomly within the three-dimensional cube [0.0,0.0,0.0]--[1.0,1.0,1.0].
Reports positions with random crystal orientations in seeds file format to STDOUT.
""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-N', dest='N', type='int', \
                  help='number of seed points to distribute [%default]')
parser.add_option('-g','--grid', dest='grid', type='int', nargs=3, \
                  help='min a,b,c grid of hexahedral box %default')
parser.add_option('-r', '--rnd', dest='randomSeed', type='int', \
                  help='seed of random number generator [%default]')

parser.set_defaults(randomSeed = 0)
parser.set_defaults(grid = [16,16,16])
parser.set_defaults(N = 20)

(options, extras) = parser.parse_args()

Npoints = reduce(lambda x, y: x * y, options.grid)
if 0 in options.grid: 
  file['croak'].write('invalid grid a b c.\n')
  sys.exit()
if options.N > Npoints: 
  sys.stderr.write('Warning: more seeds than grid points at minimum resolution.\n')
  options.N = Npoints

seeds = numpy.zeros((3,options.N),float)
numpy.random.seed(options.randomSeed)

grainEuler = numpy.random.rand(3,options.N)
grainEuler[0,:] *= 360.0
grainEuler[1,:] = numpy.arccos(2*grainEuler[1,:]-1)*180.0/math.pi
grainEuler[2,:] *= 360.0

seedpoint = numpy.random.permutation(Npoints)[:options.N]
seeds[0,:]=(numpy.mod(seedpoint                                   ,options.grid[0])\
                                                            +numpy.random.random())/options.grid[0]
seeds[1,:]=(numpy.mod(seedpoint//                  options.grid[0],options.grid[1])\
                                                            +numpy.random.random())/options.grid[1]
seeds[2,:]=(numpy.mod(seedpoint//(options.grid[1]*options.grid[0]),options.grid[2])\
                                                            +numpy.random.random())/options.grid[2]

print "5\theader"
print "$Id$"
print "grid\ta %i\tb %i\tc %i"%(options.grid[0],options.grid[1],options.grid[2],)
print "microstructures\t%i"%options.N
print "randomSeed\t%i"%(options.randomSeed)
print "x\ty\tz\tphi1\tPhi\tphi2"

numpy.savetxt(sys.stdout,numpy.transpose(numpy.concatenate((seeds,grainEuler),axis = 0)),fmt='%10.6f',delimiter='\t')
