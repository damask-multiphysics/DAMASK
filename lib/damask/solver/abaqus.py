# -*- coding: UTF-8 no BOM -*-

# $Id$

from .solver import Solver


class Abaqus(Solver):


  def __init__(self,version='',solver=''):                                                          # example version string: 6.12-2, solver: std or exp
    self.solver='Abaqus'
    if version =='':
      import subprocess
      process = subprocess.Popen(['abaqus', 'information=release'],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
      self.version = process.stdout.readlines()[1].split()[1]
      print self.version
    else:
      self.version = version
    
    if solver.lower() in ['','std','standard']:
      self.solver = 'std'
    elif solver.lower() in ['exp','explicit']:
      self.solver = 'exp'
    else:
      raise Exception('unknown Abaqus solver %'%solver)

  def return_run_command(self,model):
    import subprocess
    import re
    import damask
    env=damask.Environment()
    shortVersion = re.sub('[\.,-]', '',self.version)
    try:
      cmd='abq'+shortVersion
      subprocess.check_output(['abq'+shortVersion,'information=release'])
    except OSError:                                                                                 # link to abqXXX not existing
      cmd='abaqus'
      process = subprocess.Popen(['abaqus','information=release'],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
      detectedVersion = process.stdout.readlines()[1].split()[1]
      if self.version != detectedVersion:
        raise Exception('found Abaqus version %s, but requested %s'%(detectedVersion,self.version))
    return '%s -job %s -user %s/code/DAMASK_abaqus_%s interactive'%(cmd,model,env.rootDir(),self.solver)
      


