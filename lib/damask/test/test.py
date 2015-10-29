# -*- coding: UTF-8 no BOM -*-

# $Id$

import os, sys, shlex, inspect
import subprocess,shutil,string
import logging, logging.config
import damask
from optparse import OptionParser

class Test():
  '''
     General class for testing.
     Is sub-classed by the individual tests.
  '''

  variants = []
  
  def __init__(self,test_description):

    logger = logging.getLogger()
    logger.setLevel(0)
    fh = logging.FileHandler('test.log')                                                            # create file handler which logs even debug messages
    fh.setLevel(logging.DEBUG)
    full = logging.Formatter('%(asctime)s - %(levelname)s: \n%(message)s')
    fh.setFormatter(full)
    ch = logging.StreamHandler(stream=sys.stdout)                                                   # create console handler with a higher log level
    ch.setLevel(logging.INFO)
# create formatter and add it to the handlers
    plain = logging.Formatter('%(message)s')
    ch.setFormatter(plain)
# add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    logging.info('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n' \
         +'----------------------------------------------------------------\n' \
         +'| '+test_description+'\n' \
         +'----------------------------------------------------------------')
    self.dirBase = os.path.dirname(os.path.realpath(sys.modules[self.__class__.__module__].__file__))
    self.parser = OptionParser(
    description = test_description+' (using class: $Id run_test.py 1285 2012-02-09 08:54:09Z MPIE\m.diehl $)',
    usage='./test.py [options]')
    self.updateRequested = False

  def execute(self):
    '''
    Run all variants and report first failure.
    '''
    if not self.testPossible(): return -1
    self.clean()
    self.prepareAll()
    for variant in xrange(len(self.variants)):
      try:
        self.prepare(variant)
        self.run(variant)
        self.postprocess(variant)
        if self.updateRequested:                                                                         # update requested
          self.update(variant)
        elif not self.compare(variant):                                                                  # no update, do comparison
          return variant+1                                                                               # return culprit
      except Exception as e :
        logging.critical('\nWARNING:\n %s\n'%e)
        return variant+1                                                                                 # return culprit
    return 0
  
  def testPossible(self):
    '''
    Check if test is possible or not (e.g. no license available).
    '''
    return True
    
  def clean(self):
    '''
    Delete directory tree containing current results.
    '''
    status = True

    try:
      shutil.rmtree(self.dirCurrent())
    except:
      logging.warning('removal of directory "%s" not possible...'%(self.dirCurrent()))
      status = status and False

    try:
      os.mkdir(self.dirCurrent())
    except:
      logging.critical('creation of directory "%s" failed...'%(self.dirCurrent()))
      status = status and False

    return status
    
  def prepareAll(self):
    '''
    Do all necessary preparations for the whole test
    '''
    return True

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
    logging.debug('Update not necessary')
    return True


  def dirReference(self):
    '''
    Directory containing reference results of the test.
    '''
    return os.path.normpath(os.path.join(self.dirBase,'reference/'))


  def dirCurrent(self):
    '''
    Directory containing current results of the test.
    '''
    return os.path.normpath(os.path.join(self.dirBase,'current/'))

  
  def dirProof(self):
    '''
    Directory containing human readable proof of correctness for the test.
    '''
    return os.path.normpath(os.path.join(self.dirBase,'proof/'))

    
  def fileInReference(self,file):
    '''
    Path to a file in the refrence directory for the test.
    '''
    return os.path.join(self.dirReference(),file)


  def fileInCurrent(self,file):
    '''
    Path to a file in the current results directory for the test.
    '''
    return os.path.join(self.dirCurrent(),file)

  
  def fileInProof(self,file):
    '''
    Path to a file in the proof directory for the test.
    '''
    return os.path.join(self.dirProof(),file)

    
  def copy_Reference2Current(self,sourcefiles=[],targetfiles=[]):
    
    if len(targetfiles) == 0: targetfiles = sourcefiles
    for i,file in enumerate(sourcefiles):
      try:
        shutil.copy2(self.fileInReference(file),self.fileInCurrent(targetfiles[i]))  
      except:
        logging.critical('Reference2Current: Unable to copy file %s'%file)

 
  def copy_Base2Current(self,sourceDir,sourcefiles=[],targetfiles=[]):
    
    source=os.path.normpath(os.path.join(self.dirBase,'../../'+sourceDir))
    if len(targetfiles) == 0: targetfiles = sourcefiles
    for i,file in enumerate(sourcefiles):
      try:
        shutil.copy2(os.path.join(source,file),self.fileInCurrent(targetfiles[i]))  
      except:
        logging.error(os.path.join(source,file))
        logging.critical('Base2Current: Unable to copy file %s'%file)


  def copy_Current2Reference(self,sourcefiles=[],targetfiles=[]):
    
    if len(targetfiles) == 0: targetfiles = sourcefiles
    for i,file in enumerate(sourcefiles):
      try:
        shutil.copy2(self.fileInCurrent(file),self.fileInReference(targetfiles[i]))  
      except:
        logging.critical('Current2Reference: Unable to copy file %s'%file)

        
  def copy_Proof2Current(self,sourcefiles=[],targetfiles=[]):
    
    if len(targetfiles) == 0: targetfiles = sourcefiles
    for i,file in enumerate(sourcefiles):
      try:
        shutil.copy2(self.fileInProof(file),self.fileInCurrent(targetfiles[i]))  
      except:
        logging.critical('Proof2Current: Unable to copy file %s'%file)

        
  def copy_Current2Current(self,sourcefiles=[],targetfiles=[]):
    
    for i,file in enumerate(sourcefiles):
      try:
        shutil.copy2(self.fileInReference(file),self.fileInCurrent(targetfiles[i]))  
      except:
        logging.critical('Current2Current: Unable to copy file %s'%file)


  def execute_inCurrentDir(self,cmd,streamIn=None):

    logging.info(cmd)
    out,error = damask.util.execute(cmd,streamIn,self.dirCurrent())

    logging.info(error)
    logging.debug(out)
    
    return out,error
    

    
  def compare_Array(self,File1,File2):

    import numpy as np
    logging.info('comparing\n '+File1+'\n '+File2)
    table1 = damask.ASCIItable(File1,readonly=True)
    table1.head_read()
    len1=len(table1.info)+2
    table2 = damask.ASCIItable(File2,readonly=True)
    table2.head_read()
    len2=len(table2.info)+2

    refArray = np.nan_to_num(np.genfromtxt(File1,missing_values='n/a',skip_header = len1,autostrip=True))
    curArray = np.nan_to_num(np.genfromtxt(File2,missing_values='n/a',skip_header = len2,autostrip=True))

    if len(curArray) == len(refArray):
      refArrayNonZero = refArray[refArray.nonzero()]
      curArray = curArray[refArray.nonzero()]
      max_err=np.max(abs(refArrayNonZero[curArray.nonzero()]/curArray[curArray.nonzero()]-1.))
      max_loc=np.argmax(abs(refArrayNonZero[curArray.nonzero()]/curArray[curArray.nonzero()]-1.))
      refArrayNonZero = refArrayNonZero[curArray.nonzero()]
      curArray = curArray[curArray.nonzero()]
      print(' ********\n * maximum relative error %e for %e and %e\n ********'
                  %(max_err, refArrayNonZero[max_loc],curArray[max_loc]))
      return max_err
    else:
       raise Exception('mismatch in array size to compare')


  def compare_ArrayRefCur(self,ref,cur=''):
    
    if cur =='': cur = ref
    refName = self.fileInReference(ref)
    curName = self.fileInCurrent(cur)
    return self.compare_Array(refName,curName)

    
  def compare_ArrayCurCur(self,cur0,cur1):
    
    cur0Name = self.fileInCurrent(cur0)
    cur1Name = self.fileInCurrent(cur1)
    return self.compare_Array(cur0Name,cur1Name)

  def compare_Table(self,headings0,file0,headings1,file1,normHeadings='',normType=None,\
                                     absoluteTolerance=False,perLine=False,skipLines=[]):
    
    import numpy as np
    logging.info('comparing ASCII Tables\n %s \n %s'%(file0,file1))
    if normHeadings == '': normHeadings = headings0

    if len(headings0) == len(headings1) == len(normHeadings):                                         #check if comparison is possible and determine lenght of columns
      dataLength = len(headings0)
      length       = [1   for i in xrange(dataLength)]
      shape        = [[]  for i in xrange(dataLength)]
      data         = [[]  for i in xrange(dataLength)]
      maxError     = [0.0 for i in xrange(dataLength)]
      absTol       = [absoluteTolerance for i in xrange(dataLength)]
      column       = [[1 for i in xrange(dataLength)] for j in xrange(2)]
 
      norm         = [[]  for i in xrange(dataLength)]
      normLength   = [1   for i in xrange(dataLength)]
      normShape    = [[]  for i in xrange(dataLength)]
      normColumn   = [1   for i in xrange(dataLength)]

      for i in xrange(dataLength):
        if headings0[i]['shape'] != headings1[i]['shape']: 
          raise Exception('shape mismatch when comparing %s with %s '%(headings0[i]['label'],headings1[i]['label']))
        shape[i] = headings0[i]['shape']
        for j in xrange(np.shape(shape[i])[0]):
          length[i] *= shape[i][j]
        normShape[i] = normHeadings[i]['shape']
        for j in xrange(np.shape(normShape[i])[0]):
          normLength[i] *= normShape[i][j]
    else:
      raise Exception('trying to compare %i with %i normed by %i data sets'%(len(headings0),len(headings1),len(normHeadings)))

    table0 = damask.ASCIItable(file0,readonly=True)
    table0.head_read()
    table1 = damask.ASCIItable(file1,readonly=True)
    table1.head_read()   

    for i in xrange(dataLength):
      key0 =    {True :'1_%s',
                 False:'%s'   }[length[i]>1]%headings0[i]['label']
      key1 =    {True :'1_%s',
                 False:'%s'   }[length[i]>1]%headings1[i]['label']
      normKey = {True :'1_%s',
                 False:'%s'   }[normLength[i]>1]%normHeadings[i]['label']
      if key0 not in table0.labels:
        raise Exception('column %s not found in 1. table...\n'%key0)
      elif key1 not in table1.labels:
        raise Exception('column %s not found in 2. table...\n'%key1)
      elif normKey not in table0.labels:
        raise Exception('column %s not found in 1. table...\n'%normKey)
      else:
        column[0][i]  = table0.labels.index(key0)                                  # remember columns of requested data
        column[1][i]  = table1.labels.index(key1)                                  # remember columns of requested data in second column
        normColumn[i] = table0.labels.index(normKey)                               # remember columns of requested data in second column
    
    line0 = 0
    while table0.data_read():                                                     # read next data line of ASCII table
      if line0 not in skipLines:
        for i in xrange(dataLength):
          myData = np.array(map(float,table0.data[column[0][i]:\
                                                     column[0][i]+length[i]]),'d')
          normData = np.array(map(float,table0.data[normColumn[i]:\
                                                       normColumn[i]+normLength[i]]),'d')
          data[i] = np.append(data[i],np.reshape(myData,shape[i]))
          if normType == 'pInf':
            norm[i] = np.append(norm[i],np.max(np.abs(normData)))
          else:
            norm[i] = np.append(norm[i],np.linalg.norm(np.reshape(normData,normShape[i]),normType))
      line0 +=1
    
    for i in xrange(dataLength):
      if not perLine: norm[i] = [np.max(norm[i]) for j in xrange(line0-len(skipLines))]
      data[i] = np.reshape(data[i],[line0-len(skipLines),length[i]])
      if any(norm[i]) == 0.0 or absTol[i]:
        norm[i] = [1.0 for j in xrange(line0-len(skipLines))]
        absTol[i] = True
        if perLine:
          logging.warning('At least one norm of %s in 1. table is 0.0, using absolute tolerance'%headings0[i]['label'])
        else:
          logging.warning('Maximum norm of %s in 1. table is 0.0, using absolute tolerance'%headings0[i]['label'])

    line1 = 0
    while table1.data_read():                                                     # read next data line of ASCII table
      if line1 not in skipLines:
        for i in xrange(dataLength):
          myData = np.array(map(float,table1.data[column[1][i]:\
                                                     column[1][i]+length[i]]),'d')
          maxError[i] = max(maxError[i],np.linalg.norm(np.reshape(myData-data[i][line1-len(skipLines),:],shape[i]))/
                                                                                   norm[i][line1-len(skipLines)])
      line1 +=1

    if (line0 != line1): raise Exception('found %s lines in 1. table and %s in 2. table'%(line0,line1))

    logging.info(' ********')
    for i in xrange(dataLength):
      if absTol[i]:
        logging.info(' * maximum absolute error %e for %s and %s'%(maxError[i],headings0[i]['label'],headings1[i]['label']))
      else:
        logging.info(' * maximum relative error %e for %s and %s'%(maxError[i],headings0[i]['label'],headings1[i]['label']))
    logging.info(' ********')
    return maxError


  def compare_Table2(self,file0,file1,headings0=None,headings1=None,rtol=1e-5,atol=1e-8,threshold = -1.0,debug=False):
    
    '''
      compare two tables with np.allclose
      threshold can be used to ignore small values (put any negative number to disable)
      table will be row-wise normalized
    '''
    #http://stackoverflow.com/questions/8904694/how-to-normalize-a-2-dimensional-numpy-array-in-python-less-verbose
    import numpy as np
    logging.info('comparing ASCII Tables\n %s \n %s'%(file0,file1))

    if headings1 == None: headings1=headings0

    table0 = damask.ASCIItable(file0,readonly=True)
    table0.head_read()
    table0.data_readArray(headings0)
    row_sums0 = table0.data.sum(axis=1)*table0.data.shape[0]
    table0.data /= row_sums0[:,np.newaxis]

    table1 = damask.ASCIItable(file1,readonly=True)
    table1.head_read()   
    table1.data_readArray(headings1)
    row_sums1 = table1.data.sum(axis=1)*table1.data.shape[0]
    table1.data /= row_sums1[:,np.newaxis]
    
    if debug:
      t0 = np.where(np.abs(table0.data)<threshold,0.0,table0.data)
      t1 = np.where(np.abs(table1.data)<threshold,0.0,table1.data)
      print np.amin(np.abs(t1)*rtol+atol-np.abs(t0-t1))
      i = np.argmin(np.abs(t1)*rtol+atol-np.abs(t0-t1))
      print t0.flatten()[i],t1.flatten()[i]

    return np.allclose(np.where(np.abs(table0.data)<threshold,0.0,table0.data),
                       np.where(np.abs(table1.data)<threshold,0.0,table1.data),rtol,atol)


    
  def compare_TableRefCur(self,headingsRef,ref,headingsCur='',cur='',normHeadings='',normType=None,\
                                                 absoluteTolerance=False,perLine=False,skipLines=[]):
    
    if cur == '': cur = ref
    if headingsCur == '': headingsCur = headingsRef
    refName = self.fileInReference(ref)
    curName = self.fileInCurrent(cur)
    return self.compare_Table(headingsRef,refName,headingsCur,curName,normHeadings,normType,
                                                               absoluteTolerance,perLine,skipLines)

    
  def compare_TableCurCur(self,headingsCur0,Cur0,Cur1,headingsCur1='',normHeadings='',normType=None,\
                                                 absoluteTolerance=False,perLine=False,skipLines=[]):
    
    if headingsCur1 == '': headingsCur1 = headingsCur0
    cur0Name = self.fileInCurrent(Cur0)
    cur1Name = self.fileInCurrent(Cur1)
    return self.compare_Table(headingsCur0,cur0Name,headingsCur1,cur1Name,normHeadings,normType,
                                                               absoluteTolerance,perLine,skipLines)


  def report_Success(self,culprit):
    
    if culprit == 0:
      logging.critical('%s passed.'%({False: 'The test',
                         True: 'All %i tests'%(len(self.variants))}[len(self.variants) > 1]))
      logging.critical('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
      return 0
    if culprit == -1:
      logging.warning('Warning: Could not start test')
      return 0
    else:
      logging.critical(' ********\n * Test %i failed...\n ********'%(culprit))
      logging.critical('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
      return culprit
