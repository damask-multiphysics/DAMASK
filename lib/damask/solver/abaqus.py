# -*- coding: UTF-8 no BOM -*-

# $Id$

from .solver import Solver


class Abaqus(Solver):

  def __init__(self,version='',method=''):
    self.solver='Abaqus'
    if version =='':
      self.version = version
      cmd = "abaqus information=release"
      import subprocess
      process = subprocess.Popen(cmd,stdout = subprocess.PIPE,stderr = subprocess.PIPE,shell=True)
      self.version = process.stdout.readlines()[1].split()[1]
    else:
      self.version = version
    
    if method.lower() in ['','std','standard']:
      self.method = 'std'
    elif method.lower() in ['exp','explicit']:
      self.method = 'exp'
    else:
      self.method = 'chose either std or exp'

  def return_run_command(self,model):
    import subprocess
    import re
    import damask
    env=damask.Environment()
    shortVersion = re.sub('[\.,-]', '',self.version)
    try:
      subprocess.Popen(['abq'+shortVersion,'information=release'])
      cmd='abq'+shortVersion
    except subprocess.CalledProcessError:
      process = subprocess.Popen(cmd,stdout = subprocess.PIPE,stderr = subprocess.PIPE,shell=True)
      cmd='abaqus'
      if self.version != process.stdout.readlines()[1].split()[1]: raise Exception
    return '%s -job %s -user %s/code/DAMASK_abaqus_%s interactive'%(cmd,model,env.rootDir(),self.method)
      


