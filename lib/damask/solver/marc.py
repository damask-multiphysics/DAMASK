# $Id$

from .solver import Solver


class Marc(Solver):

#--------------------------
  def __init__(self):
#--------------------------
    self.solver = 'Marc'
    self.releases = { \
              '2010.2':['linux64',''],
              '2010':  ['linux64',''],
              '2008r1':[''],
              '2007r1':[''],
              '2005r3':[''],
             }


#--------------------------
  def version(self,rootRelation = ''):
#--------------------------
    import os,damask.environment

    MSCpath = damask.environment.Environment(rootRelation).pathInfo['msc']
    
    for release,subdirs in sorted(self.releases.items(),reverse=True):
      for subdir in subdirs:
        libPath = '%s/mentat%s/shlib/%s'%(MSCpath,release,subdir)
        if os.path.exists(libPath): return release
        else: continue
    
    return ''
      
    
#--------------------------
  def libraryPath(self,rootRelation = ''):
#--------------------------
    import os,damask.environment

    MSCpath = damask.environment.Environment(rootRelation).pathInfo['msc']
    
    for release,subdirs in sorted(self.releases.items(),reverse=True):
      for subdir in subdirs:
        libPath = '%s/mentat%s/shlib/%s'%(MSCpath,release,subdir)
        if os.path.exists(libPath): return libPath
        else: continue
    
    return ''
  

#--------------------------
  def submit_job(self,
#--------------------------
                 rootRelation = '',
                 run_marc_path='/msc/marc2010/tools/',
                 subroutine_dir=None,
                 subroutine_name='DAMASK_marc2010',
                 compile=False,
                 compiled_dir='../../../code/',
                 modelname='one_element_model',
                 jobname='job1',
                 #IOdir='',
                 host=[]
                 ):
    import os,damask.environment
    import subprocess, shlex   
    import shutil
    
    damaskEnv = damask.environment.Environment(rootRelation)
    if subroutine_dir  is None: subroutine_dir  = damaskEnv.binDir()
    if subroutine_name is None: subroutine_name = 'DAMASK_marc' + self.version(rootRelation)
    if run_marc_path   is None: run_marc_path   = os.path.join(damaskEnv.pathInfo['msc'],self.version(rootRelation),'tools/')

    # Define all options [see Marc Installation and Operation Guide, pp 23]
    run_marc = os.path.jion(run_marc_path,'run_marc')
    jid = ' -jid ' + modelname + '_' + jobname
    compilation=' -u ' + subroutine_dir + subroutine_name + '.f90'+' -save y'
    options=' -nprocd 1  -autorst 0 -ci n  -cr n  -dcoup 0 -b no -v no'
    cmd=run_marc+jid+options
          
    if compile:
      cmd += compilation
      print 'job submission with compilation.'
    else:
      shutil.copy2(subroutine_dir+subroutine_name+'.f90','./'+subroutine_name+'.f90')
      shutil.copy2(compiled_dir+subroutine_name+'.marc','./'+subroutine_name+'.marc')
      prog = ' -prog ' + subroutine_name
      cmd += prog
      print 'Job submission without compilation, using %s'%prog
    out = open('out.log','w')
    print(cmd)
    #print shlex.split(cmd)
    self.p = subprocess.Popen(shlex.split(cmd),stdout=out,stderr=subprocess.STDOUT)
    self.p.wait()
    out.close()
      
#--------------------------
  def exit_number_from_outFile(self,outFile=None):
#--------------------------
    import string
    fid_out = open(outFile,'r')
    for ln in fid_out:
      if (string.find(ln,'tress iteration') is not -1):
        print ln
      elif (string.find(ln,'Exit number') is not -1):
        substr = ln[string.find(ln,'Exit number'):len(ln)]
        exitnumber = substr[12:16]
        fid_out.close()
        return int(exitnumber)
    fid_out.close()            
    return -1    
