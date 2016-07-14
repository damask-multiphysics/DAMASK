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
          items = l.split('=')
          if len(items) == 2: 
            self.options[items[0].upper()] = \
              re.sub('\$\{*DAMASK_ROOT\}*',self.rootDir(),os.path.expandvars(items[1]))             # expand all shell variables and DAMASK_ROOT
      
  def isAvailable(self,software,Nneeded =-1):
    licensesNeeded = {'abaqus'  :5,
                      'standard':5,
                      'explicit':5}
    if Nneeded == -1: Nneeded = licensesNeeded[software]
    try:
      cmd = """ ssh mulicense2 "/Stat_Flexlm | grep 'Users of %s: ' | cut -d' ' -f7,13" """%software
      process = subprocess.Popen(shlex.split(cmd),stdout = subprocess.PIPE,stderr = subprocess.PIPE)
      licenses = map(int, process.stdout.readline().split())
      try:
        if licenses[0]-licenses[1] >= Nneeded:
          return 0
        else:
          print('%s missing licenses for %s'%(licenses[1] + Nneeded - licenses[0],software))
          return licenses[1] + Nneeded - licenses[0]
      except IndexError:
        print('Could not retrieve license information for %s'%software)
        return 127
    except:
      return 126
