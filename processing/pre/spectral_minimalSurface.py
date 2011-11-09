#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,math
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


minimal_surfaces = ['gyroid','d']

surface = {
            'gyroid': lambda x,y,z: math.sin(x)*math.cos(y)+math.sin(y)*math.cos(z)+math.cos(x)*math.sin(z),
            'd':      lambda x,y,z: math.cos(z)*math.sin(x+y)+math.sin(z)*math.cos(x-y),
          }
  

parser = OptionParser(option_class=extendedOption, usage='%prog', description = """
Produce VTK point file from geom data
""" + string.replace('$Id: spectral_geomCheck 994 2011-09-05 13:38:10Z MPIE\p.eisenlohr $','\n','\\n')
)

parser.add_option('-t', '--type', dest='type', type='string', \
                  help='type of minimal surface [%s]'%(','.join(minimal_surfaces)))
parser.add_option('-c', '--threshold', dest='threshold', type='float', \
                  help='threshold value defining minimal surface')
parser.add_option('-p', '--periods', dest='periods', type='int', \
                  help='number of repetitions of unit cell')
parser.add_option('-r', '--resolution', dest='resolution', type='int', nargs=3, \
                  help='resolution in x, y, and z')
parser.add_option('-d', '--dimension', dest='dimension', type='float', nargs=3, \
                  help='physical dimension in x, y, and z')

parser.set_defaults(type = 'gyroid')
parser.set_defaults(threshold = 0.5)
parser.set_defaults(periods = 1)
parser.set_defaults(resolution = [8,8,8])
parser.set_defaults(dimension = [1.0,1.0,1.0])

(options, args) = parser.parse_args()

print "3 header"
print "resolution\ta %i\tb %i\tc %i"%(options.resolution[0],options.resolution[1],options.resolution[2],)
print "dimension\tx %g\ty %g\tz %g"%(options.dimension[0],options.dimension[1],options.dimension[2],)
print "homogenization  1"

for z in range(options.resolution[2]):
  for y in range(options.resolution[1]):
    for x in range(options.resolution[0]):
      print {True:'1',False:'2'}\
            [options.threshold > 
             surface[options.type](options.periods*2.0*math.pi*x/options.resolution[0],
                                   options.periods*2.0*math.pi*y/options.resolution[1],
                                   options.periods*2.0*math.pi*z/options.resolution[2],
                                  )]
