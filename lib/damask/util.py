# -*- coding: UTF-8 no BOM -*-

# damask utility functions
import threading
import numpy as np
from optparse import OptionParser, Option

# -----------------------------
# Matlab like trigonometric functions that take and return angles in degrees.
# -----------------------------
for f in ['cos', 'sin', 'tan']:
  exec('def %sd(deg): return (np.%s(np.deg2rad(deg)))'%(f,f))
  exec('def a%sd(val): return (np.rad2deg(np.arc%s(val)))'%(f,f))


# -----------------------------
def gridLocation(idx,res):
# -----------------------------
  return ( idx  % res[0], \
         ( idx // res[0]) % res[1], \
         ( idx // res[0] // res[1]) % res[2] )


# -----------------------------
def gridIndex(location,res):
# -----------------------------
  return ( location[0] % res[0]                   + \
         ( location[1] % res[1]) * res[0]          + \
         ( location[2] % res[2]) * res[1] * res[0]   )


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


# -----------------------------
class backgroundMessage(threading.Thread):
# -----------------------------
    
    def __init__(self):
        threading.Thread.__init__(self)
        self.message = ''
        self.new_message = ''
        self.counter = 0
        self.symbols = ['- ', '\ ', '| ', '/ ']
        self.waittime = 0.5
    
    def __quit__(self):
        length = len(self.message) + len(self.symbols[self.counter])
        sys.stderr.write(chr(8)*length + ' '*length + chr(8)*length)
        sys.stderr.write('')
    
    def run(self):
        while not threading.enumerate()[0]._Thread__stopped:
            time.sleep(self.waittime)
            self.update_message()
        self.__quit__()

    def set_message(self, new_message):
        self.new_message = new_message
        self.print_message()
    
    def print_message(self):
        length = len(self.message) + len(self.symbols[self.counter])
        sys.stderr.write(chr(8)*length + ' '*length + chr(8)*length)                                # delete former message
        sys.stderr.write(self.symbols[self.counter] + self.new_message)                             # print new message
        self.message = self.new_message
        
    def update_message(self):
        self.counter = (self.counter + 1)%len(self.symbols)
        self.print_message()

