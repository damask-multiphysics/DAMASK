#!/usr/bin/env python
# $Id$

import os, sys
import subprocess,shutil
import damask

class Test():
  '''
     General class for testing.
     Is sub-classed by the individual tests.
  '''

  variants = []
  
  def __init__(self):
    self.dirBase = os.path.dirname(os.path.realpath(sys.argv[0]))
  
  def execute(self,variants = [],update = []):
    '''
    Run all variants and report first failure.
    '''
    if len(variants) == 0: variants = xrange(len(self.variants))       # iterate over all variants

    self.clean()
    for variant in variants:
      self.prepare(variant)
      self.run(variant)
      self.postprocess(variant)
      if variant in update:
        self.update(variant)
      elif not self.compare(variant):
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
      print('removal of directory "%s" failed...'%(self.dirCurrent()))
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
      shutil.copy2(self.fileInReference(file),self.fileInCurrent(file))  

  def copy_Current2Current(self,files=[]):
    for file in files:
      shutil.copy2(self.fileInCurrent(file[0]),self.fileInCurrent(file[1]))  

