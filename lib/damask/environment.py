# $Id$

import os,sys,string,re,subprocess,shlex

class Environment():
  __slots__ = [ \
                'rootRelation',
                'pathInfo',
              ]

  def __init__(self,rootRelation = '.'):
    self.rootRelation = rootRelation                    
    self.options = {}
    self.get_options()

  def relPath(self,relative = '.'):
    return os.path.join(self.rootDir(),relative)

  def rootDir(self):
    damask_root = os.getenv('DAMASK_ROOT')
    if damask_root == '' or damask_root == None:            # env not set
      if sys.argv[0] == '':                                 # no caller path
        cwd = os.getcwd()
      else:
        cwd = {False: os.path.dirname(sys.argv[0]),
                True:                 sys.argv[0]}[os.path.isdir(sys.argv[0])]
      damask_root = os.path.normpath(os.path.join(os.path.realpath(cwd),self.rootRelation))

    return damask_root

  def binDir(self):
    damask_bin  = os.getenv('DAMASK_BIN')
    if damask_bin  == '' or damask_bin  == None: damask_bin = self.relPath('bin/')
    return damask_bin

  def get_options(self):
    try:                                                    # check for user-defined pathinfo
      file = open(self.relPath('installation/options'))
      content = map(lambda string: string.strip(),file.readlines())
      file.close()
      for line in content:
        if not (line.startswith('#') or line == ''):
          items = line.split('=') + ['','']
          if items[1] != '':                                # nothing specified
            self.options[items[0].upper()] = items[1]
    except:
      pass
      
  def isAvailable(self,software,noNeeded=-1):
    licensesNeeded = {'abaqus'  :5,
                      'standard':5,
                      'explicit':5}
    if noNeeded == -1: noNeeded = licensesNeeded[software]
    cmd = """ ssh mulicense2 "/Stat_Flexlm | grep 'Users of %s: ' | cut -d' ' -f7,13" """%software
    process = subprocess.Popen(shlex.split(cmd),stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    licenses = map(int, process.stdout.readline().split())
    try:
      if licenses[0]-licenses[1] >= noNeeded:
        return 0
      else:
        print(licenses[1] + noNeeded - licenses[0], 'missing licenses for %s'%software)
        return licenses[1] + noNeeded - licenses[0]
    except IndexError:
      print('Could not retrieve license information for %s'%software)
      return 127
