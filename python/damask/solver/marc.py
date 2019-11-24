import os
import subprocess
import shlex

from .solver import Solver
import damask

class Marc(Solver):
  """Wrapper to run DAMASK with MSCMarc."""

  def __init__(self,version=damask.Environment().options['MARC_VERSION']):
    """
    Create a Marc solver object.
    
    Parameters
    ----------
    version : float
        Marc version

    """
    self.solver  ='Marc'
    self.version = damask.environment.Environment().options['MARC_VERSION']


#--------------------------
  def libraryPath(self):

    path_MSC = damask.environment.Environment().options['MSC_ROOT']   
    path_lib = '{}/mentat{}/shlib/linux64'.format(path_MSC,self.version)

    return path_lib if os.path.exists(path_lib) else ''


#--------------------------
  def toolsPath(self):

    path_MSC   = damask.environment.Environment().options['MSC_ROOT']
    path_tools = '{}/marc{}/tools'.format(path_MSC,self.version)

    return path_tools if os.path.exists(path_tools) else ''


#--------------------------
  def submit_job(self,
                 model,
                 job          = 'job1',
                 logfile      = False,
                 compile      = False,
                 optimization ='',
                ):


    damaskEnv = damask.environment.Environment()
   
    user = os.path.join(damaskEnv.relPath('src'),'DAMASK_marc{}.{}'.format(self.version,'f90' if compile else 'marc'))
    if not os.path.isfile(user):
      raise FileNotFoundError("DAMASK4Marc ({}) '{}' not found".format(('source' if compile else 'binary'),user))

    # Define options [see Marc Installation and Operation Guide, pp 23]
    script = 'run_damask_{}mp'.format(optimization)
    
    cmd = os.path.join(self.toolsPath(),script) + \
          ' -jid ' + model + '_' + job + \
          ' -nprocd 1  -autorst 0 -ci n  -cr n  -dcoup 0 -b no -v no'

    if compile: cmd += ' -u ' + user + ' -save y'
    else:       cmd += ' -prog ' + os.path.splitext(user)[0]

    print('job submission with{} compilation: {}'.format('' if compile else 'out',user))
    if logfile: log = open(logfile, 'w')
    print(cmd)
    process = subprocess.Popen(shlex.split(cmd),stdout = log,stderr = subprocess.STDOUT)
    log.close()
    process.wait()
      
#--------------------------
  def exit_number_from_outFile(self,outFile=None):
    import string
    exitnumber = -1
    fid_out = open(outFile,'r')
    for line in fid_out:
      if (string.find(line,'tress iteration') != -1):
        print(line)
      elif (string.find(line,'Exit number')   != -1):
        substr = line[string.find(line,'Exit number'):len(line)]
        exitnumber = int(substr[12:16])

    fid_out.close()
    return exitnumber
