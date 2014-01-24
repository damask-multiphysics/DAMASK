# $Id$

import os,sys,string,re,subprocess,shlex

class Environment():
  __slots__ = [ \
                'rootRelation',
                'options',
              ]

  def __init__(self,rootRelation = '.'):
    self.rootRelation = rootRelation
    self.options = {}
    self.get_options()

  def relPath(self,relative = '.'):
    return os.path.join(self.rootDir(),relative)

  def rootDir(self):
    return os.path.normpath(os.path.join(os.path.realpath(__file__),'../../../'))

  def get_options(self):
    if os.path.isfile(os.path.join(os.getenv('HOME'),'.damask/damask.conf')):
      configFile = os.path.join(os.getenv('HOME'),'.damask/damask.conf')
    else:
      configFile = '/etc/damask.conf'
    with open(self.relPath(configFile)) as file:
      for line in file:
        l = line.strip()
        if not (l.startswith('#') or l == ''):
          items = l.split('=') + ['','']
          if items[1] != '':                                # nothing specified
            self.options[items[0].upper()] = items[1]
      
  def isAvailable(self,software,Nneeded =-1):
    licensesNeeded = {'abaqus'  :5,
                      'standard':5,
                      'explicit':5}
    if Nneeded == -1: Nneeded = licensesNeeded[software]
    cmd = """ ssh mulicense2 "/Stat_Flexlm | grep 'Users of %s: ' | cut -d' ' -f7,13" """%software
    process = subprocess.Popen(shlex.split(cmd),stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    licenses = map(int, process.stdout.readline().split())
    try:
      if licenses[0]-licenses[1] >= Nneeded:
        return 0
      else:
        print(licenses[1] + Nneeded - licenses[0], 'missing licenses for %s'%software)
        return licenses[1] + Nneeded - licenses[0]
    except IndexError:
      print('Could not retrieve license information for %s'%software)
      return 127
