#!/usr/bin/env python
import os, sys
import subprocess,shutil

import damask_tools
import msc_tools

damask_tools.DAMASK_TOOLS().check_env()
  
class DAMASK_TEST():
    '''
       General class for testing.
       Is sub-classed by the individual tests.
    '''
    modelname=None
    jobname=None
    testdir=None
    spectral_options=None
    compile=False
    
    has_variants=False  # False ==> A single test case is run
    #has_variants=True   # True ==> Need to define the test_variants generator

    def test_variants(self):
        ''' 
           If has_subtests == True this method defines the generator for subtests
           This generator must be defined in each test, 
           depending on what to change: orientations, parameters,....
           Below is an EXAMPLE.
        '''
        subtest_orientations=[[0.0,90.0,0.0],[0.0,0.0,90.0]]
        for i,o in enumerate(subtest_orientations):
          from damask_tools import MATERIAL_CONFIG
          mat=MATERIAL_CONFIG()
          mat.read('material.config_base')
          mat.add_texture(label='Grain001',
                          type ='(gauss)',
                          eulers = o)
          mat.write(overwrite=True)
          print(mat.data['texture']['Grain001'])
          testlabel='orientation_%03i'%i
          yield(testlabel)

    def run_test(self):
        res=[]
        if self.has_variants:
          for t in self.subtests():
            print(t)
            val=self.run_single_test()
            res.append(val==True)
        else:
          val=self.run_single_test()
          res.append(val==True)
        if all(res) is True:
          return True     
        print(res)
        return False 
          
    def run_single_test(self):
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
        
    def calc_current_results(self):
        '''
           Should be defined in the individual tests
        '''
        pass

    def calc_marc(self,compile=None):
        '''
           Calculates current results for MSC.Marc
        '''
        if compile is None: compile=self.compile
        self.copy_from_ref=[self.modelname+'_'+self.jobname+'.dat',
                   self.modelname+'.mfd', # for dev
                   'material.config'
        ]
        self.copy_files_from_reference_results()
        os.chdir('./current_results')                
        m=msc_tools.MSC_TOOLS()
        m.submit_job(compile=compile, compiled_dir='../../../code/')
        print('simulation submitted')
        self.exit_number=m.exit_number_from_outFile(outFile=self.modelname+'_'+self.jobname+'.out')
       
        if not self.exit_number==3004:
          print('Job did not exit with No. 3004')
          return False 
        return True
        
    def calc_spectral(self, compile=None):
        pass

    def copy_files_from_reference_results(self):
        for file in self.copy_from_ref:
            shutil.copy2('./reference_results/%s'%file,'./current_results/%s'%file)  
            # Note: possibly symlinking? No, because copy is OS independent.
        
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
    
