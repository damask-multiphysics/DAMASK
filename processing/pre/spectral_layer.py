#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,math,numpy,random
from optparse import OptionParser, OptionGroup, Option


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

parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
Construct continuous Microstructure (e.g. from EBSD) with Layer
""" + string.replace('$Id: spectral_randomSeeding.py 1423 2012-03-31 12:42:49Z MPIE\m.diehl $','\n','\\n')
)

parser.add_option('-r','--res', dest='res', type='int', nargs=2, \
                  help='Continuous Fourier Points in x, y [%default]')
parser.add_option('-l','--layer', dest='layer', type='int', nargs=2, \
                  help='Constant Layer of Fourier Points x, y [%default]')

parser.set_defaults(res=[16,16])
parser.set_defaults(layer=[0,0])

(options, filenames) = parser.parse_args()

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    files.append({'name':name,'output':open(name+'_tmp','w')})
for file in files:
  if file['name'] != 'STDIN': print file['name']
  file['output'].write("{0:1d} {1:6s}\n".format(1,'header'))
  file['output'].write("{0:s} {1:8d} {2:s} {3:8d} {4:s}\n".format('resolution a',options.res[0]+options.layer[0],'b',options.res[1]+options.layer[1],'c 1'))
  for y in xrange(options.res[1]+options.layer[1]):
    file['output'].write('%i copies of %i\n'%(options.layer[0],options.res[0]*options.res[1]+2))
    if (y<options.layer[1]): 
      file['output'].write('%i copies of %i\n'%(options.res[0],options.res[0]*options.res[1]+1))
    else:
      file['output'].write('%i to %i\n'%((y-options.layer[1])*options.res[0]+1,(y-options.layer[1])*options.res[0]+options.res[0]))
  file['output'].close()
