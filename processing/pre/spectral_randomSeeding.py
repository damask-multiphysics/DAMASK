#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math,numpy,random
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

identifiers = {
        'resolution': ['a','b','c'],
        'dimension':  ['x','y','z'],
          }
mappings = {
        'resolution': lambda x: int(x),
        'dimension':  lambda x: float(x),
          }


parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
Offset microstructure index for points which see a microstructure different from themselves within a given (cubic) vicinity,
i.e. within the region close to a grain/phase boundary.
""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-f', '--file', dest='filename', type="string", \
                  help='output seed file name [%default]')
parser.add_option('-s', '--seed', dest='randomSeed', type='int', \
                  help='seed of random number generator [%default]')
parser.add_option('-n', '--ngrains', dest='N_Seeds', type='int', \
                  help='seed of random number generator[%default]')
parser.add_option('-r','--res', dest='res', type='int', nargs=3, \
                  help='Min Fourier points in x, y, z [%default]')

parser.set_defaults(filename = 'seeds')
parser.set_defaults(randomSeed = 0)
parser.set_defaults(res=[16,16,16])
parser.set_defaults(N_Seeds=50)

(options, filenames) = parser.parse_args()

if options.N_Seeds > options.res[0]*options.res[1]*options.res[2]: 
  print 'Warning: Number of grains exceeds min resolution'
  options.N_Seeds = options.res[0]*options.res[1]*options.res[2]

seeds = numpy.zeros((3,options.N_Seeds),float)
numpy.random.seed(options.randomSeed)

grainEuler=numpy.random.rand(3,options.N_Seeds)
grainEuler[0,:] = 360*grainEuler[0,:]
grainEuler[1,:] = numpy.arccos(2*grainEuler[1,:]-1)*180/math.pi
grainEuler[2,:] = 360*grainEuler[2,:]

seedpoint = numpy.random.permutation(options.res[0]*options.res[1]*options.res[2])[:options.N_Seeds]
seeds[0,:]=(numpy.mod(seedpoint                                 ,options.res[0])+numpy.random.random())/options.res[0]
seeds[1,:]=(numpy.mod(seedpoint//                options.res[0] ,options.res[1])+numpy.random.random())/options.res[1]
seeds[2,:]=(numpy.mod(seedpoint//(options.res[1]*options.res[0]),options.res[2])+numpy.random.random())/options.res[2]

f = open(options.filename+'.seeds', 'w')

f.write("{0:1d} {1:6s}\n".format(4,'header'))
f.write("{0:s} {1:8d} {2:s} {3:8d} {4:s} {5:8d}\n".format('resolution a',options.res[0],'b',options.res[1],'c',options.res[2]))
f.write("{0:s} {1:8d}\n".format('grains',options.N_Seeds))
f.write("{0:s} {1:8d}\n".format('random seed',options.randomSeed))
f.write("x y z phi1 Phi phi2\n")
f.close()
f=file(options.filename+'.seeds','a')
numpy.savetxt(f,numpy.transpose(numpy.concatenate((seeds,grainEuler),axis=0)),fmt='%10.6f',delimiter=' ')
f.close()
