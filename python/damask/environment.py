# -*- coding: UTF-8 no BOM -*-

import os,subprocess,shlex,re

class Environment():
  __slots__ = [ \
                'options',
              ]

  def __init__(self):
    self.options = {}
    self.get_options()

  def relPath(self,relative = '.'):
    return os.path.join(self.rootDir(),relative)

  def rootDir(self):
    return os.path.normpath(os.path.join(os.path.realpath(__file__),'../../../'))

  def get_options(self):
    with open(self.relPath(self.rootDir()+'/CONFIG')) as configFile:
      for line in configFile:
        l = re.sub('^set ', '', line).strip()                                                       # remove "set" (tcsh) when setting variables
        if l and not l.startswith('#'):
          items = re.split(r'\s*=\s*',l)
          if len(items) == 2: 
            self.options[items[0].upper()] = \
              re.sub('\$\{*DAMASK_ROOT\}*',self.rootDir(),os.path.expandvars(items[1]))             # expand all shell variables and DAMASK_ROOT
