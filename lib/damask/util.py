# -*- coding: UTF-8 no BOM -*-

# damask utility functions
import math
from optparse import OptionParser, Option

try:
  import numpy
except:
  numpy=None

# Matlab like trigonometric functions that take and return angles in degrees.
for f in ['cos', 'sin', 'tan']:
    if numpy:
      exec('def %sd(deg): return (numpy.%s(numpy.deg2rad(deg)))'%(f,f))
      exec('def a%sd(val): return (numpy.rad2deg(numpy.arc%s(val)))'%(f,f))
    else:  
      exec('def %sd(deg): return (math.%s(deg/180.*math.pi))'%(f,f))
      exec('def a%sd(val): return (math.a%s(val)*180./math.pi)'%(f,f))

# -----------------------------
class extendableOption(Option):
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

