#!/usr/bin/env python
import os, sys
import subprocess,shutil

import damask_tools; reload(damask_tools)
import msc_tools; reload(msc_tools)  

damask_tools.DAMASK_TOOLS().check_env()
  
class DAMASK_TEST():
    modelname=None
    jobname=None
    testdir=None
    spectral_options=None
    orientations=[]
    compile=False
        
    def run_test(self):
        self.modelname='one_element_model'    
        self.jobname='job1'
        self.testdir='2001_hex_plastic'
        self.orientations=[]

        self.clean_current_results()        
        if self.calc_current_results() is False:
          return False
        print('simulation finished')          
        self.postprocess()
        if self.compare_to_reference() is False:
          return False
        print 'Test OK'
        return True

    def clean_current_results(self):
        try:
          shutil.rmtree('current_results')
        except:
          print('Could not delete current_results')
        os.mkdir('current_results')

    def calc_current_results(self,compile=None):
        #theDir = os.path.split(sys.argv[0])[0]
        #os.chdir(theDir)
        #os.chdir('..')
        #os.chdir('%s/testing'%os.getenv('DAMASK_ROOT'))
        if compile is None: compile=self.compile
        self.copy_from_ref=[self.modelname+'_'+self.jobname+'.dat',
                   self.modelname+'.mfd', # for dev
                   'material.config'
        ]
        for file in self.copy_from_ref:
            shutil.copy2('./reference_results/%s'%file,'./current_results/%s'%file)  
            # Note: possibly symlinking? No, because copy is OS independent.

        os.chdir('./current_results')

        m=msc_tools.MSC_TOOLS()
        m.submit_job(compile=compile, compiled_dir='../../../code/')
        print('simulation submitted')
        self.exit_number=m.exit_number_from_outFile(outFile=self.modelname+'_'+self.jobname+'.out')
       
        if not self.exit_number==3004:
          print('Job did not exit with No. 3004')
          return False 
        return True  
        
    def postprocess(self):
        #print 'postprocessing results ...'
        #os.system('%s/processing/post/postResults --es "Comp 33 of Stress" %s.t16 --range 100 100 1'%(os.getenv('DAMASK_ROOT'),self.modelname+'_'+self.jobname))
        print 'postprocessing results ...'
        file=open('../postprocessing.cmd','r')
        postproc=file.readlines()
        file.close()
        for cmd in postproc:    # PHILIP: suggestion to just execute the script "postprocessing" directly within a shell, i.e. os.system('../postprocessing')
          print(cmd)
          os.system(cmd)        # PHILIP: reason is that for loops and the like get broken with line by line execution from here...
# CLAUDIO: Actually that's what we had before - I stole the code from one of your scripts because for lengthy postprocessing, the user can then see the progress. I don't get the part about the breaking loops, let's discuss tomorrow. 


    def compare_to_reference(self,tol=1e-5):    
        import string
        print 'comparing results against reference_results...'  
        txt_file=self.modelname+'_'+self.jobname+'.txt'
        cur=self.read_val_from_file(fname='postProc/'+txt_file)
        ref=self.read_val_from_file(fname='../reference_results/postProc/'+txt_file)
        
        err=abs((ref/cur)-1.) # relative tolerance
        #err=abs(ref-cur)      # absolute tolerance
        
        if err>tol: 
          print 'Current value:   %e'%cur
          print 'Reference value: %e'%ref
          print('err: %e > tol: %e'%(err,tol))
          return False
        print('err: %e < tol: %e'%(err,tol))  
        return True
        
    def read_val_from_file(self,fname=None):
        fid=open(fname,'r')
        rl=fid.readlines()
        print rl
        cur=rl[-1]
        lst=cur.split('\t')
        print lst
        val=float(lst[-1].rstrip())
        print val
        fid.close()
        return val
        

if __name__ == "__main__":
    test=DAMASK_TESTER()
    test.run_test()
    
