import subprocess

from .solver import Solver
import damask


class Abaqus(Solver):
  """Wrapper to run DAMASK with Abaqus."""

  def __init__(self,version=int(damask.Environment().options['ABAQUS_VERSION'])):
    """
    Create a Abaqus solver object.
    
    Parameters
    ----------
    version : integer
        Abaqus version

    """
    self.solver  ='Abaqus'
    self.version = int(version)

  def return_run_command(self,model):
    env=damask.Environment()
    try:
      cmd='abq'+self.version
      subprocess.check_output([cmd,'information=release'])
    except OSError:                                                                                 # link to abqXXX not existing
      cmd='abaqus'
      process = subprocess.Popen(['abaqus','information=release'],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
      detectedVersion = int(process.stdout.readlines()[1].split()[1].decode('utf-8'))
      if self.version != detectedVersion:
        raise Exception('found Abaqus version {}, but requested {}'.format(detectedVersion,self.version))
    return '{} -job {} -user {}/src/DAMASK_abaqus interactive'.format(cmd,model,env.rootDir())
