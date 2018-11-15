# -*- coding: UTF-8 no BOM -*-


from .solver import Solver


class Marc(Solver):

  def __init__(self):
    self.solver = 'Marc'


#--------------------------
  def version(self):
    import damask.environment

    return damask.environment.Environment().options['MARC_VERSION']
      
    
#--------------------------
  def libraryPath(self,release = ''):
    import os,damask.environment

    MSCpath     = damask.environment.Environment().options['MSC_ROOT']
    if len(release) == 0: release = self.version()
    
    path = '{}/mentat{}/shlib/linux64'.format(MSCpath,release)

    return path if os.path.exists(path) else ''


#--------------------------
  def toolsPath(self,release = ''):
    import os,damask.environment

    MSCpath = damask.environment.Environment().options['MSC_ROOT']
    if len(release) == 0: release = self.version()

    path = '%s/marc%s/tools'%(MSCpath,release)

    return path if os.path.exists(path) else ''
  

#--------------------------
  def submit_job(self,
                 release      = '',
                 model        = 'model',
                 job          = 'job1',
                 logfile      = None,
                 compile      = False,
                 optimization ='',
                ):

    import os,damask.environment
    import subprocess,shlex
    
    if len(release) == 0: release = self.version()


    damaskEnv = damask.environment.Environment()
    
    user = os.path.join(damaskEnv.relPath('src/'),'DAMASK_marc')                                    # first guess, might be updated in the following lines
    if compile:
      if os.path.isfile(os.path.join(damaskEnv.relPath('src/'),'DAMASK_marc%s.f90'%release)):
        user = os.path.join(damaskEnv.relPath('src/'),'DAMASK_marc%s'%release)
    else:
      if os.path.isfile(os.path.join(damaskEnv.relPath('src/'),'DAMASK_marc%s.marc'%release)):
        user = os.path.join(damaskEnv.relPath('src/'),'DAMASK_marc%s'%release)

    # Define options [see Marc Installation and Operation Guide, pp 23]
    script = 'run_damask_%smp'%({False:'',True:optimization}[optimization!=''])
    
    cmd = os.path.join(self.toolsPath(release),script) + \
          ' -jid ' + model + '_' + job + \
          ' -nprocd 1  -autorst 0 -ci n  -cr n  -dcoup 0 -b no -v no'

    if compile: cmd += ' -u ' + user+'.f90' + ' -save y'
    else:       cmd += ' -prog ' + user

    print('job submission with%s compilation: %s'%({False:'out',True:''}[compile],user))
    if logfile:
      log = open(logfile, 'w')
    print(cmd)
    self.p = subprocess.Popen(shlex.split(cmd),stdout = log,stderr = subprocess.STDOUT)
    log.close()
    self.p.wait()
      
#--------------------------
  def exit_number_from_outFile(self,outFile=None):
    import string
    exitnumber = -1
    fid_out = open(outFile,'r')
    for ln in fid_out:
      if (string.find(ln,'tress iteration') is not -1):
        print(ln)
      elif (string.find(ln,'Exit number') is not -1):
        substr = ln[string.find(ln,'Exit number'):len(ln)]
        exitnumber = int(substr[12:16])

    fid_out.close()
    return exitnumber
