# -*- coding: UTF-8 no BOM -*-

from .solver import Solver
import damask
import subprocess

class Abaqus(Solver):

  def __init__(self,version=''):                                                                    # example version string: 2017
    self.solver='Abaqus'
    if version =='':
      version = damask.Environment().options['ABAQUS_VERSION']
    else:
      self.version = version

  def return_run_command(self,model):
    env=damask.Environment()
    try:
      cmd='abq'+self.version
      subprocess.check_output([cmd,'information=release'])
    except OSError:                                                                                 # link to abqXXX not existing
      cmd='abaqus'
      process = subprocess.Popen(['abaqus','information=release'],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
      detectedVersion = process.stdout.readlines()[1].split()[1].decode('utf-8')
      if self.version != detectedVersion:
        raise Exception('found Abaqus version {}, but requested {}'.format(detectedVersion,self.version))
    return '{} -job {} -user {}/src/DAMASK_abaqus interactive'.format(cmd,model,env.rootDir())
