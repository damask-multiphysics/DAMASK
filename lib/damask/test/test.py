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
    
  def compare_TableRefCur(self,headings,ref,cur=''):
    
    if cur =='': cur = ref
    if len(headings) == 1 or type(headings) == dict: headings += headings
    if len(headings) != 2 : print 'headings should be array of length 1 or 2 containing dictionaries or a dictionary'
    if len(headings[0]) != len(headings[1]): print 'mismatch in headings to compare'
    
    refName = self.fileInReference(ref)
    curName = self.fileInCurrent(cur)
    return self.compare_Table(headings,refName,curName)
    
  def compare_Table(self,headings,File1,File2):

    print 'comparing ASCII Tables\n' , File1,'\n', File2
    table1 = damask.ASCIItable(open(File1))
    table1.head_read()
    table2 = damask.ASCIItable(open(File2))
    table2.head_read()   
    datainfo = {                                                               # list of requested labels per datatype
                'scalar':     {'len':1,
                               'label1':[],
                               'label2':[]},
                'vector':     {'len':3,
                               'label1':[],
                               'label2':[]},
                'tensor':     {'len':9,
                               'label1':[],
                               'label2':[]},
                }

    for label in headings[0]:
      datainfo[headings[0][label]]['label1'] += [label]
    for label in headings[1]:
      datainfo[headings[1][label]]['label2'] += [label]

    active = [{},{}]
    column = [{},{}]

    for datatype,info in datainfo.items():
      for label in info['label1']:
        key = {True :'1_%s',
               False:'%s'   }[info['len']>1]%label
        if key not in table1.labels:
          sys.stderr.write('column %s not found in 1. table...\n'%key)
        else:
          if datatype not in active[0]: active[0][datatype] = []
          if datatype not in column[0]: column[0][datatype] = {}
          active[0][datatype].append(label)
          column[0][datatype][label] = table1.labels.index(key)                   # remember columns of requested data
      for label in info['label2']:
        key = {True :'1_%s',
               False:'%s'   }[info['len']>1]%label
        if key not in table2.labels:
          sys.stderr.write('column %s not found in 2. table...\n'%key)
        else:
          if datatype not in active[1]: active[1][datatype] = []
          if datatype not in column[1]: column[1][datatype] = {}
          active[1][datatype].append(label)
          column[1][datatype][label] = table2.labels.index(key)                   # remember columns of requested data
    
    while table1.data_read():                                                     # read next data line of ASCII table
      for datatype,labels in active[0].items():                                   # loop over vector,tensor
        for label in labels:                                                      # loop over all requested norms
          data1[label] += numpy.array(map(float,table1.data[column[0][datatype][label]:
                              column[0][datatype][label]+datainfo[datatype]['len']]),'d').reshape(3,3)
    print myData
  def report_Success(self,culprit):
    
    if culprit < 0:
      print '%s passed.'%({False: 'The test',
                         True: 'All %i tests'%(len(self.variants))}[len(self.variants) > 1])
    else:
     print ' ********\n * Test %i failed...\n ********'%(culprit+1)

    print '\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
