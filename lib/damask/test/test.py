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
    
    print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n' \
         +'----------------------------------------------------------------\n' \
         +'| '+test_description+'\n' \
         +'----------------------------------------------------------------'
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
    if len(update)   == 0 and self.options.update: print ' This test has no reference to update'
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

  def copy_Reference2Current(self,sourcefiles=[],targetfiles=[]):
    
    if len(targetfiles) == 0: targetfiles = sourcefiles
    for i,file in enumerate(sourcefiles):
      try:
        shutil.copy2(self.fileInReference(file),self.fileInCurrent(targetfiles[i]))  
      except:
        print 'Reference2Current: Unable to copy file ', file
 
  def copy_Current2Reference(self,sourcefiles=[],targetfiles=[]):
    
    if len(targetfiles) == 0: targetfiles = sourcefiles
    for i,file in enumerate(sourcefiles):
      try:
        shutil.copy2(self.fileInCurrent(file),self.fileInReference(targetfiles[i]))  
      except:
        print 'Current2Reference: Unable to copy file ', file
 
  def copy_Current2Current(self,sourcefiles=[],targetfiles=[]):
    
    for i,file in enumerate(sourcefiles):
      try:
        shutil.copy2(self.fileInReference(file),self.fileInCurrent(targetfiles[i]))  
      except:
        print 'Current2Current: Unable to copy file ', file

  def execute_inCurrentDir(self,cmd,outfile='execute_log.txt'):
    
    os.chdir(self.dirCurrent())
    file=open(outfile,'a+')
    print cmd
    process = subprocess.Popen(shlex.split(cmd),stdout = file,stderr = subprocess.STDOUT)
    process.wait()
    file.close()
    
  def compare_Array(self,File1,File2):
  
    import numpy
    print 'comparing\n ' , File1,'\n ', File2
    refFile = open(File1)
    table = damask.ASCIItable(refFile)
    table.head_read()
    refFile.close()
    refArray = numpy.nan_to_num(numpy.genfromtxt(File1,missing_values='n/a',skip_header = len(table.info)+1,autostrip=True))
    curArray = numpy.nan_to_num(numpy.genfromtxt(File2,missing_values='n/a',skip_header = len(table.info)+1,autostrip=True))
    if len(curArray) ==  len(refArray):
      refArrayNonZero = refArray[refArray.nonzero()]
      curArray = curArray[refArray.nonzero()]
      max_err=numpy.max(abs(refArrayNonZero[curArray.nonzero()]/curArray[curArray.nonzero()]-1.))
      max_loc=numpy.argmax(abs(refArrayNonZero[curArray.nonzero()]/curArray[curArray.nonzero()]-1.))
      refArrayNonZero = refArrayNonZero[curArray.nonzero()]
      curArray = curArray[curArray.nonzero()]
      print ' ********\n * maximum relative error ',max_err,' for ', refArrayNonZero[max_loc],' and ',curArray[max_loc],'\n ********'
      return max_err
    else:
      print ' ********\n * mismatch in array size to compare \n ********'
      return sys.float_info.max
      
  def compare_ArrayRefCur(self,ref,cur=''):
    
    if cur =='': cur = ref
    refName = self.fileInReference(ref)
    curName = self.fileInCurrent(cur)
    return self.compare_Array(refName,curName)
    
  def compare_TableRefCur(self,headingsRef,ref,headingsCur='',cur=''):
    
    if headingsCur == '': headingsCur = headingsRef
    if cur == '': cur = ref
    refName = self.fileInReference(ref)
    curName = self.fileInCurrent(cur)
    return self.compare_Table(headingsRef,refName,headingsCur,curName)
    
  def compare_Table(self,headings0,file0,headings1,file1):
    
    import numpy
    print 'comparing ASCII Tables\n' , file0,'\n', file1

    if len(headings0) == len(headings1):                                                             #check if comparison is possible and determine lenght of columns
      dataLength = len(headings0)
      length       = [1   for i in xrange(dataLength)]
      shape        = [[]  for i in xrange(dataLength)]
      data         = [[]  for i in xrange(dataLength)]
      maxError     = [0.0 for i in xrange(dataLength)]
      maxNorm      = [0.0 for i in xrange(dataLength)]

      column = [[length],[length]]
      for i in xrange(dataLength):
        if headings0[i]['shape'] != headings1[i]['shape']: 
          raise Exception('shape mismatch when comparing ', headings0[i]['label'], ' with ', headings1[i]['label'])
        shape[i] = headings0[i]['shape'] 
        for j in xrange(numpy.shape(headings0[i]['shape'])[0]):
          length[i] = headings0[i]['shape'][j]
    else:
      raise Exception('trying to compare ', len(headings0), ' with ', len(headings1), ' data sets')


    table0 = damask.ASCIItable(open(file0))
    table0.head_read()
    table1 = damask.ASCIItable(open(file1))
    table1.head_read()   

    for i in xrange(dataLength):
      key0 = {True :'1_%s',
             False:'%s'   }[length[i]>1]%headings0[i]['label']
      key1 = {True :'1_%s',
             False:'%s'   }[length[i]>1]%headings1[i]['label']
       
      if key0 not in table0.labels:
        raise Exception('column %s not found in 1. table...\n'%key0)
      elif key1 not in table1.labels:
        raise Exception('column %s not found in 2. table...\n'%key1)
      else:
        column[0][i] = table0.labels.index(key0)                   # remember columns of requested data
        column[1][i] = table1.labels.index(key1)                   # remember columns of requested data in second column

    line0 = 0
    while table0.data_read():                                                     # read next data line of ASCII table
      line0 +=1
      for i in xrange(dataLength):
        myData = numpy.array(map(float,table0.data[column[0][i]:\
                                                   column[0][i]+length[i]]),'d')
        maxNorm[i] = max(maxNorm[i],numpy.linalg.norm(numpy.reshape(myData,shape[i])))
        data[i]=numpy.append(data[i],myData)
    
    for i in xrange(dataLength):
      data[i] = numpy.reshape(data[i],[line0,length[i]])
    
    line1 = 0
    while table1.data_read():                                                     # read next data line of ASCII table
      for i in xrange(dataLength):
        myData = numpy.array(map(float,table1.data[column[1][i]:\
                                                   column[1][i]+length[i]]),'d')
        maxError[i] = max(maxError[i],numpy.linalg.norm(numpy.reshape(myData-data[i][line1,:],shape[i])))
        line1 +=1
        
    if (line0 != line1): raise Exception('found ', line0, ' lines in 1. table and ', line1, ' in 2. table')
    
    for i in xrange(dataLength):
      maxError[i] = maxError[i]/maxNorm[i]
    
    return maxError
      
  def report_Success(self,culprit):
    
    if culprit < 0:
      print '%s passed.'%({False: 'The test',
                         True: 'All %i tests'%(len(self.variants))}[len(self.variants) > 1])
    else:
     print ' ********\n * Test %i failed...\n ********'%(culprit+1)

    print '\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
