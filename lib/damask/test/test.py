#!/usr/bin/env python
import os, sys
import subprocess,shutil

import damask
import msc_tools

class Test():
    '''
       General class for testing.
       Is sub-classed by the individual tests.
    '''
    #define those according to your test
    modelname=None
    jobname=None
    test_dir=None
    spectral_options=None
    compile=False = None
    post_txt = None
    tol = 0.0
    
    has_variants=False  # False ==> A single test case is run
    #has_variants=True   # True ==> Need to define the test_variants generator

    def test_variants(self):
        ''' 
           If has_subtests == True this method defines the generator for test variants
           This generator must be defined in each test, 
           depending on what to change: orientations, parameters,....
           Below is an EXAMPLE.
        '''
        #maybe get rid of has_variants through testing for None return value...
        yield(None)

    def test_variants_example(self):
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
          for t in self.test_variants():
            print '###############################################'
            print '###############################################'
            print(t)
            print '###############################################'
            val=self.run_single_test(t)
            res.append(val==True)
        else:
          val=self.run_single_test()
          res.append(val==True)
        if all(res) is True:
          return True     
        print(res)
        return False 
          
    def run_single_test(self,variant):
        self.clean_current_results()        
        if self.calc_current_results(variant) is False:
          return False
        print('simulation finished')          
        self.postprocess()
        if self.compare_to_reference() is False:
          print '++++++++ Test not OK +++++++++'
          return False
        print 'Test OK'
        return True

    def clean_current_results(self):
        os.chdir(self.test_dir)
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
        self.copy_from_ref_list()
        self.copy_files_from_reference_results()
        os.chdir('./current_results')                
        #m=msc_tools.MSC_TOOLS()
        #m.submit_job(compile=compile, compiled_dir='../../../code/')
        damask.solver.Marc().submit_job(compile=compile, compiled_dir='../../../code/')
        print('simulation submitted')
        self.exit_number=m.exit_number_from_outFile(outFile=self.modelname+'_'+self.jobname+'.out')
       
        if not self.exit_number==3004:
          print('Job did not exit with No. 3004')
          return False 
        return True
        
    def calc_spectral(self, compile=None):
        pass

    def copy_from_ref_list(self):
        self.copy_from_ref=[self.modelname+'_'+self.jobname+'.dat',
                            self.modelname+'.mfd', # for dev
                            'material.config',
                           ]

        
    def copy_files_from_reference_results(self):
        for file in self.copy_from_ref:
            shutil.copy2('./reference_results/%s'%file,'./current_results/%s'%file)  
            # Note: possibly symlinking? No, because copy is OS independent.
        
    def read_val_from_file(self,fname=None):
        file = open(fname,'r')
        headerlength = int(file.readline().split()[0]) + 1
        file.close
        import numpy as N
        val = N.loadtxt(fname,skiprows=headerlength)
        return val

    def compare_to_reference(self):    
        import string, numpy as N
        print 'comparing results against reference_results...'
        os.chdir(os.path.join(self.test_dir,'current_results'))
        cur=self.read_val_from_file(fname='postProc/'+self.post_txt)
        ref=self.read_val_from_file(fname='../reference_results/postProc/'+self.post_txt) 
        
        err=abs((ref/cur)-1.) # relative tolerance
        #err=abs(ref-cur)      # absolute tolerance
        print 'tol', self.tol
        print 'max error', N.max(err)
        if N.max(err)>self.tol:
          return False 
        return True
        
    def postprocess(self):
        print 'postprocessing results ...'
        os.chdir(self.test_dir)
        file=open('./postprocessing.cmd','r')
        postproc=file.readlines()
        file.close()
        os.chdir(os.path.join(self.test_dir,'current_results'))
        for cmd in postproc:    # PHILIP: suggestion to just execute the script "postprocessing" directly within a shell, i.e. os.system('../postprocessing')
           print(cmd)
           os.system(cmd)        # PHILIP: reason is that for loops and the like get broken with line by line execution from here...

    
