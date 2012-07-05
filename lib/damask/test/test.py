#!/usr/bin/env python
# $Id$

import os, sys, shlex
import subprocess,shutil,string
import damask
from optparse import OptionParser

class Test():
  '''
     General class for testing.
     Is sub-classed by the individual tests.
  '''

  variants = []
  
  def __init__(self,test_description):
    print '\n'+test_description
    self.dirBase = os.path.dirname(os.path.realpath(sys.argv[0]))
    self.parser = OptionParser(
    description = 'Using: $Id run_test.py 1285 2012-02-09 08:54:09Z MPIE\m.diehl $',
    usage='run_test.py [options]')
  
    self.parser.add_option("-u", "--update", action="store_true",\
                                    dest="update",\
                                    help="use current test results as new reference")
    self.parser.set_defaults(update = False)
  
    (self.options, self.args) = self.parser.parse_args()
    
  def execute(self,variants = [],update = []):
    
    '''
    Run all variants and report first failure.
    '''
    if len(variants) == 0: variants = xrange(len(self.variants))       # iterate over all variants
    self.clean()
    for variant in variants:
      try:
        self.prepare(variant)
        self.run(variant)
        self.postprocess(variant)
        if variant in update:
          self.update(variant)
        elif not self.compare(variant):
          return variant
      except Exception,e :
        print '\nWARNING:\n %s\n'%e
        return variant
    return -1

  def clean(self):
    '''
    Delete directory tree containing current results.
    '''
    status = True

    try:
      shutil.rmtree(self.dirCurrent())
    except:
      print('removal of directory "%s" not possible...'%(self.dirCurrent()))
      status = status and False

    try:
      os.mkdir(self.dirCurrent())
    except:
      print('creation of directory "%s" failed...'%(self.dirCurrent()))
      status = status and False

    return status


  def prepare(self,variant):
    '''
    Do all necessary preparations for the run of each test variant
    '''
    return True
  

  def run(self,variant):
    '''
    Execute the requested test variant.
    '''
    return True


  def postprocess(self,variant):
    '''
    Perform post-processing of generated results for this test variant.
    '''
    return True


  def compare(self,variant):
    '''
    Compare reference to current results.
    '''
    return True


  def update(self,variant):
    '''
    Update reference with current results.
    '''
    return True


  def dirReference(self):
    '''
    '''
    return os.path.normpath(os.path.join(self.dirBase,'reference/'))

  def dirCurrent(self):
    '''
    '''
    return os.path.normpath(os.path.join(self.dirBase,'current/'))

  def fileInReference(self,file):
    '''
    '''
    return os.path.join(self.dirReference(),file)

  def fileInCurrent(self,file):
    '''
    '''
    return os.path.join(self.dirCurrent(),file)

  def copy_Reference2Current(self,files=[]):
    for file in files:
      try:
        shutil.copy2(self.fileInReference(file),self.fileInCurrent(file))  
      except:
        print 'Reference2Current: Unable to copy file ', file

  def copy_Current2Current(self,files=[]):
    for file in files:
      try:
        shutil.copy2(self.fileInCurrent(file[0]),self.fileInCurrent(file[1]))  
      except:
        print 'Current2Current: Unable to copy file ', file

  def execute_inCurrentDir(self,cmd,outfile='execute_log.txt'):
    os.chdir(self.dirCurrent())
    file=open(outfile,'a+')
    print cmd
    process = subprocess.Popen(shlex.split(cmd),stdout = file,stderr = subprocess.STDOUT)
    file.close()
    process.wait()


  def compare_Array(self,File1,File2):
    import numpy

    refFile = open(File1)
    table = damask.ASCIItable(refFile)
    table.head_read()
    refFile.close()
    refArray = numpy.nan_to_num(numpy.genfromtxt(File1,missing_values='n/a',skip_header = len(table.info)+1))
    curArray = numpy.nan_to_num(numpy.genfromtxt(File2,missing_values='n/a',skip_header = len(table.info)+1))
    print len(refArray)
    print refArray
    print len(curArray)
    print curArray
    refNonZero = refArray[refArray.nonzero()]
    curNonZero = curArray[curArray.nonzero()]
    print refNonZero
    print curNonZero
    err = abs((refNonZero/curNonZero)-1.)                                     # relative tolerance
    print err
    max_err = numpy.max(err)
    print ' ********\n * maximum relative error',max_err,'\n ********'
    return max_err
    
  def compare_ArrayRefCur(self,ref,cur):
    refName = self.fileInReference(ref)
    curName = self.fileInCurrent(cur)
    return compare_Array(refName,curName)

  def report_Success(self,culprit):
    if culprit < 0:
      print '%s passed.'%({False: 'The test',
                         True: 'All %i tests'%(len(self.variants))}[len(self.variants) > 1])
    else:
     print ' ********\n * Test %i failed...\n ********'%(culprit+1)
