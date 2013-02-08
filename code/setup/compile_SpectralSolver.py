#!/usr/bin/env python

import os,string
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP
import damask

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

    
parser = OptionParser(option_class=extendedOption, usage='%prog [options] datafile[s]', description = """
Writes version specific files for different MARC releases, adjustes the make file for the spectral solver and optionally compiles the spectral solver

""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-c', '--compile', dest='spectralCompile', action='store_true', \
                  help='compiles the spectral solver [%default]')
parser.add_option('-o', '--options', dest='makeOptions', action='extend', type='string', \
                  metavar="KEY=VALUE", \
                  help='comma-separated list of options passed to Makefile when compiling spectral code %default')
parser.set_defaults(spectralCompile = True)
parser.set_defaults(makeOptions = ['F90=ifort'])

(options, args) = parser.parse_args()

damaskEnv = damask.Environment('../../')          # script location relative to root
baseDir = damaskEnv.relPath('code/')


# compiling spectral code
if (options.spectralCompile):
  print 'base directory:', baseDir
  os.system('make --directory %s clean'%(baseDir))
  print 'make options:',' '.join(options.makeOptions)
  os.system('make --directory %s %s'%(baseDir,' '.join(options.makeOptions)))