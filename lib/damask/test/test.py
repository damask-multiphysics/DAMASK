# -*- coding: UTF-8 no BOM -*-

import os,sys,shutil
import logging,logging.config
import damask
import numpy as np
from collections import Iterable
from optparse import OptionParser

class Test():
  """
  General class for testing.

  Is sub-classed by the individual tests.
  """

  variants = []
  
  def __init__(self, **kwargs):

    defaults = {'description': '',
                'keep':          False,
                'accept':        False,
                'updateRequest': False,
                'show':          False,
                'select':        None,
                }
    for arg in defaults.keys():
      setattr(self,arg,kwargs.get(arg) if kwargs.get(arg) else defaults[arg])
    
    fh = logging.FileHandler('test.log')                                                            # create file handler which logs even debug messages
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: \n%(message)s'))

    ch = logging.StreamHandler(stream=sys.stdout)                                                   # create console handler with a higher log level
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(message)s'))

    logger = logging.getLogger()
    logger.addHandler(fh)
    logger.addHandler(ch)
    logger.setLevel(0)

    logging.info('\n'.join(['+'*40,
                            '-'*40,
                            '| '+self.description,
                            '-'*40,
                           ]))

    self.dirBase = os.path.dirname(os.path.realpath(sys.modules[self.__class__.__module__].__file__))

    self.parser = OptionParser(option_class=damask.extendableOption,
                               description = '{} (Test class version: {})'.format(self.description,damask.version),
                               usage = './test.py [options]')
    self.parser.add_option("-k", "--keep",
                           action = "store_true",
                           dest   = "keep",
                           help   = "keep current results, just run postprocessing")
    self.parser.add_option("--ok", "--accept",
                           action = "store_true",
                           dest   = "accept",
                           help   = "calculate results but always consider test as successful")
    self.parser.add_option("-l",   "--list",
                           action = "store_true",
                           dest   = "show",
                           help   = "show all test variants without actual calculation")
    self.parser.add_option("-s",   "--select",
                           dest   = "select",
                           action = 'extend', metavar = '<string LIST>',
                           help   = "run test(s) of given name only")
    self.parser.set_defaults(keep   = self.keep,
                             accept = self.accept,
                             update = self.updateRequest,
                             show   = self.show,
                             select = self.select,
                            )

    
  def variantName(self,variant):
    """Generate name of (numerical) variant."""
    return str(variant)

  def execute(self):
    """Run all variants and report first failure."""
    if not self.options.keep:
      if not self.feasible(): return -1
      self.clean()
      self.prepareAll()

    for variant,object in enumerate(self.variants):
      name = self.variantName(variant)
      if self.options.show:
        logging.critical('{}: {}'.format(variant+1,name))
      elif self.options.select is not None \
           and not (name in self.options.select or str(variant+1) in self.options.select):
        pass
      else:
        try:
          if not self.options.keep:
            self.prepare(variant)
            self.run(variant)

          self.postprocess(variant)
          
          if self.options.update:
            if self.update(variant) != 0: logging.critical('update for "{}" failed.'.format(name))
          elif not (self.options.accept or self.compare(variant)):                                    # no update, do comparison
            return variant+1                                                                          # return culprit

        except Exception as e:
          logging.critical('exception during variant execution: "{}"'.format(str(e)))
          return variant+1                                                                            # return culprit
    return 0
  
  def feasible(self):
    """Check whether test is possible or not (e.g. no license available)."""
    return True

  def clean(self):
    """Delete directory tree containing current results."""
    try:
      shutil.rmtree(self.dirCurrent())
    except:
      logging.warning('removal of directory "{}" not possible...'.format(self.dirCurrent()))

    try:
      os.mkdir(self.dirCurrent())
      return True
    except:
      logging.critical('creation of directory "{}" failed.'.format(self.dirCurrent()))
      return False
    
  def prepareAll(self):
    """Do all necessary preparations for the whole test"""
    return True

  def prepare(self,variant):
    """Do all necessary preparations for the run of each test variant"""
    return True


  def run(self,variant):
    """Execute the requested test variant."""
    return True


  def postprocess(self,variant):
    """Perform post-processing of generated results for this test variant."""
    return True


  def compare(self,variant):
    """Compare reference to current results."""
    return True


  def update(self,variant):
    """Update reference with current results."""
    logging.critical('update not supported.')
    return 1


  def dirReference(self):
    """Directory containing reference results of the test."""
    return os.path.normpath(os.path.join(self.dirBase,'reference/'))


  def dirCurrent(self):
    """Directory containing current results of the test."""
    return os.path.normpath(os.path.join(self.dirBase,'current/'))


  def dirProof(self):
    """Directory containing human readable proof of correctness for the test."""
    return os.path.normpath(os.path.join(self.dirBase,'proof/'))


  def fileInRoot(self,dir,file):
    """Path to a file in the root directory of DAMASK."""
    return os.path.join(damask.Environment().rootDir(),dir,file)


  def fileInReference(self,file):
    """Path to a file in the refrence directory for the test."""
    return os.path.join(self.dirReference(),file)


  def fileInCurrent(self,file):
    """Path to a file in the current results directory for the test."""
    return os.path.join(self.dirCurrent(),file)


  def fileInProof(self,file):
    """Path to a file in the proof directory for the test."""
    return os.path.join(self.dirProof(),file)


  def copy(self, mapA, mapB,
                 A = [], B = []):
    """
    Copy list of files from (mapped) source to target.

    mapA/B is one of self.fileInX.
    """
    if not B or len(B) == 0: B = A

    for source,target in zip(list(map(mapA,A)),list(map(mapB,B))):
      try:
        shutil.copy2(source,target)
      except:
        logging.critical('error copying {} to {}'.format(source,target))
        raise


  def copy_Reference2Current(self,sourcefiles=[],targetfiles=[]):

    if len(targetfiles) == 0: targetfiles = sourcefiles
    for i,f in enumerate(sourcefiles):
      try:
        shutil.copy2(self.fileInReference(f),self.fileInCurrent(targetfiles[i]))
      except:
        logging.critical('Reference2Current: Unable to copy file "{}"'.format(f))
        raise


  def copy_Base2Current(self,sourceDir,sourcefiles=[],targetfiles=[]):

    source=os.path.normpath(os.path.join(self.dirBase,'../../..',sourceDir))
    if len(targetfiles) == 0: targetfiles = sourcefiles
    for i,f in enumerate(sourcefiles):
      try:
        shutil.copy2(os.path.join(source,f),self.fileInCurrent(targetfiles[i]))
      except:
        logging.error(os.path.join(source,f))
        logging.critical('Base2Current: Unable to copy file "{}"'.format(f))
        raise


  def copy_Current2Reference(self,sourcefiles=[],targetfiles=[]):

    if len(targetfiles) == 0: targetfiles = sourcefiles
    for i,f in enumerate(sourcefiles):
      try:
        shutil.copy2(self.fileInCurrent(f),self.fileInReference(targetfiles[i]))
      except:
        logging.critical('Current2Reference: Unable to copy file "{}"'.format(f))
        raise


  def copy_Proof2Current(self,sourcefiles=[],targetfiles=[]):

    if len(targetfiles) == 0: targetfiles = sourcefiles
    for i,f in enumerate(sourcefiles):
      try:
        shutil.copy2(self.fileInProof(f),self.fileInCurrent(targetfiles[i]))
      except:
        logging.critical('Proof2Current: Unable to copy file "{}"'.format(f))
        raise


  def copy_Current2Current(self,sourcefiles=[],targetfiles=[]):

    for i,f in enumerate(sourcefiles):
      try:
        shutil.copy2(self.fileInReference(f),self.fileInCurrent(targetfiles[i]))
      except:
        logging.critical('Current2Current: Unable to copy file "{}"'.format(f))
        raise

  def execute_inCurrentDir(self,cmd,streamIn=None):

    logging.info(cmd)
    out,error = damask.util.execute(cmd,streamIn,self.dirCurrent())

    logging.info(error)
    logging.debug(out)

    return out,error



  def compare_Array(self,File1,File2):

    import numpy as np
    logging.info('\n '.join(['comparing',File1,File2]))
    table1 = damask.ASCIItable(name=File1,readonly=True)
    table1.head_read()
    len1=len(table1.info)+2
    table2 = damask.ASCIItable(name=File2,readonly=True)
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
      print(' ********\n * maximum relative error {} between {} and {}\n ********'.format(max_err,
                                                                                          refArrayNonZero[max_loc],
                                                                                          curArray[max_loc]))
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

  def compare_Table(self,headings0,file0,
                         headings1,file1,
                         normHeadings='',normType=None,
                         absoluteTolerance=False,perLine=False,skipLines=[]):

    import numpy as np
    logging.info('\n '.join(['comparing ASCII Tables',file0,file1]))
    if normHeadings == '': normHeadings = headings0

# check if comparison is possible and determine lenght of columns
    if len(headings0) == len(headings1) == len(normHeadings):
      dataLength = len(headings0)
      length       = [1   for i in range(dataLength)]
      shape        = [[]  for i in range(dataLength)]
      data         = [[]  for i in range(dataLength)]
      maxError     = [0.0 for i in range(dataLength)]
      absTol       = [absoluteTolerance for i in range(dataLength)]
      column       = [[1  for i in range(dataLength)] for j in range(2)]

      norm         = [[]  for i in range(dataLength)]
      normLength   = [1   for i in range(dataLength)]
      normShape    = [[]  for i in range(dataLength)]
      normColumn   = [1   for i in range(dataLength)]

      for i in range(dataLength):
        if headings0[i]['shape'] != headings1[i]['shape']:
          raise Exception('shape mismatch between {} and {} '.format(headings0[i]['label'],headings1[i]['label']))
        shape[i] = headings0[i]['shape']
        for j in range(np.shape(shape[i])[0]):
          length[i] *= shape[i][j]
        normShape[i] = normHeadings[i]['shape']
        for j in range(np.shape(normShape[i])[0]):
          normLength[i] *= normShape[i][j]
    else:
      raise Exception('trying to compare {} with {} normed by {} data sets'.format(len(headings0),
                                                                                   len(headings1),
                                                                                   len(normHeadings)))

    table0 = damask.ASCIItable(name=file0,readonly=True)
    table0.head_read()
    table1 = damask.ASCIItable(name=file1,readonly=True)
    table1.head_read()

    for i in range(dataLength):
      key0    = ('1_' if     length[i]>1 else '') +    headings0[i]['label']
      key1    = ('1_' if     length[i]>1 else '') +    headings1[i]['label']
      normKey = ('1_' if normLength[i]>1 else '') + normHeadings[i]['label']
      if key0 not in table0.labels(raw = True):
        raise Exception('column "{}" not found in first table...\n'.format(key0))
      elif key1 not in table1.labels(raw = True):
        raise Exception('column "{}" not found in second table...\n'.format(key1))
      elif normKey not in table0.labels(raw = True):
        raise Exception('column "{}" not found in first table...\n'.format(normKey))
      else:
        column[0][i]  = table0.label_index(key0)
        column[1][i]  = table1.label_index(key1)
        normColumn[i] = table0.label_index(normKey)

    line0 = 0
    while table0.data_read():                                                  # read next data line of ASCII table
      if line0 not in skipLines:
        for i in range(dataLength):
          myData = np.array(map(float,table0.data[column[0][i]:\
                                                  column[0][i]+length[i]]),'d')
          normData = np.array(map(float,table0.data[normColumn[i]:\
                                                    normColumn[i]+normLength[i]]),'d')
          data[i] = np.append(data[i],np.reshape(myData,shape[i]))
          if normType == 'pInf':
            norm[i] = np.append(norm[i],np.max(np.abs(normData)))
          else:
            norm[i] = np.append(norm[i],np.linalg.norm(np.reshape(normData,normShape[i]),normType))
      line0 += 1

    for i in range(dataLength):
      if not perLine: norm[i] = [np.max(norm[i]) for j in range(line0-len(skipLines))]
      data[i] = np.reshape(data[i],[line0-len(skipLines),length[i]])
      if any(norm[i]) == 0.0 or absTol[i]:
        norm[i] = [1.0 for j in range(line0-len(skipLines))]
        absTol[i] = True
        if perLine:
          logging.warning('At least one norm of "{}" in first table is 0.0, using absolute tolerance'.format(headings0[i]['label']))
        else:
          logging.warning('Maximum norm of "{}" in first table is 0.0, using absolute tolerance'.format(headings0[i]['label']))

    line1 = 0
    while table1.data_read():                                                  # read next data line of ASCII table
      if line1 not in skipLines:
        for i in range(dataLength):
          myData = np.array(map(float,table1.data[column[1][i]:\
                                                  column[1][i]+length[i]]),'d')
          maxError[i] = max(maxError[i],np.linalg.norm(np.reshape(myData-data[i][line1-len(skipLines),:],shape[i]))/
                                                                                   norm[i][line1-len(skipLines)])
      line1 +=1

    if (line0 != line1): raise Exception('found {} lines in first table but {} in second table'.format(line0,line1))

    logging.info(' ********')
    for i in range(dataLength):
      if absTol[i]:
        logging.info(' * maximum absolute error {} between {} and {}'.format(maxError[i],
                                                                             headings0[i]['label'],
                                                                             headings1[i]['label']))
      else:
        logging.info(' * maximum relative error {} between {} and {}'.format(maxError[i],
                                                                             headings0[i]['label'],
                                                                             headings1[i]['label']))
    logging.info(' ********')
    return maxError


  def compare_TablesStatistically(self,
                     files = [None,None],                                      # list of file names
                     columns = [None],                                         # list of list of column labels (per file)
                     meanTol = 1.0e-4,
                     stdTol = 1.0e-6,
                     preFilter = 1.0e-9):
    """
    Calculate statistics of tables

    threshold can be used to ignore small values (a negative number disables this feature)
    """
    if not (isinstance(files, Iterable) and not isinstance(files, str)):       # check whether list of files is requested
      files = [str(files)]

    tables = [damask.ASCIItable(name = filename,readonly = True) for filename in files]
    for table in tables:
      table.head_read()

    columns += [columns[0]]*(len(files)-len(columns))                          # extend to same length as files
    columns = columns[:len(files)]                                             # truncate to same length as files

    for i,column in enumerate(columns):
      if column is None: columns[i] = tables[i].labels(raw = True)             # if no column is given, read all

    logging.info('comparing ASCIItables statistically')
    for i in range(len(columns)):
      columns[i] = columns[0]  if not columns[i] else \
                 ([columns[i]] if not (isinstance(columns[i], Iterable) and not isinstance(columns[i], str)) else \
                   columns[i]
                 )
      logging.info(files[i]+':'+','.join(columns[i]))

    if len(files) < 2: return True                                             # single table is always close to itself...

    data = []
    for table,labels in zip(tables,columns):
      table.data_readArray(labels)
      data.append(table.data)
      table.close()


    for i in range(1,len(data)):
      delta = data[i]-data[i-1]
      normBy = (np.abs(data[i]) + np.abs(data[i-1]))*0.5
      normedDelta = np.where(normBy>preFilter,delta/normBy,0.0)
      mean = np.amax(np.abs(np.mean(normedDelta,0)))
      std = np.amax(np.std(normedDelta,0))
      logging.info('mean: {:f}'.format(mean))
      logging.info('std:  {:f}'.format(std))

    return (mean<meanTol) & (std < stdTol)



  def compare_Tables(self,
                     files   = [None,None],                                    # list of file names
                     columns = [None],                                         # list of list of column labels (per file)
                     rtol    = 1e-5,
                     atol    = 1e-8,
                     debug   = False):
    """Compare multiple tables with np.allclose"""
    if not (isinstance(files, Iterable) and not isinstance(files, str)):       # check whether list of files is requested
      files = [str(files)]

    if len(files) < 2: return True                                             # single table is always close to itself...

    tables = [damask.ASCIItable(name = filename,readonly = True) for filename in files]
    for table in tables:
      table.head_read()

    columns += [columns[0]]*(len(files)-len(columns))                          # extend to same length as files
    columns = columns[:len(files)]                                             # truncate to same length as files

    for i,column in enumerate(columns):
      if column is None: columns[i] = tables[i].labels(raw = False)            # if no column is given, use all

    logging.info('comparing ASCIItables')
    for i in range(len(columns)):
      columns[i] = columns[0]  if not columns[i] else \
                 ([columns[i]] if not (isinstance(columns[i], Iterable) and not isinstance(columns[i], str)) else \
                   columns[i]
                 )
      logging.info(files[i]+': '+','.join(columns[i]))

    dimensions = tables[0].label_dimension(columns[0])                         # width of each requested column
    maximum = np.zeros_like(columns[0],dtype=float)                            # one magnitude per column entry
    data = []                                                                  # list of feature table extracted from each file (ASCII table)

    for i,(table,labels) in enumerate(zip(tables,columns)):
      if np.any(dimensions != table.label_dimension(labels)):                  # check data object consistency
        logging.critical('Table {} differs in data layout.'.format(files[i]))
        return False
      table.data_readArray(labels)                                             # read data, ...
      data.append(table.data)                                                  # ... store, ...
      table.close()                                                            # ... close

      for j,label in enumerate(labels):                                        # iterate over object labels
        maximum[j] = np.maximum(
                       maximum[j],
                       np.amax(np.linalg.norm(table.data[:,table.label_indexrange(label)],
                                              axis=1))
                      )                                                        # find maximum Euclidean norm across rows

    maximum = np.where(maximum > 0.0, maximum, 1.0)                            # avoid div by zero for zero columns
    maximum = np.repeat(maximum,dimensions)                                    # spread maximum over columns of each object

    for i in range(len(data)):
      data[i] /= maximum                                                       # normalize each table
      logging.info('shape of data {}: {}'.format(i,data[i].shape))

    if debug:
      violators = np.absolute(data[0]-data[1]) > atol + rtol*np.absolute(data[1])
      logging.info('shape of violators: {}'.format(violators.shape))
      for j,culprits in enumerate(violators):
        goodguys = np.logical_not(culprits)
        if culprits.any():
          logging.info('{} has {}'.format(j,np.sum(culprits)))
          logging.info('deviation: {}'.format(np.absolute(data[0][j]-data[1][j])[culprits]))
          logging.info('data     : {}'.format(np.absolute(data[1][j])[culprits]))
          logging.info('deviation: {}'.format(np.absolute(data[0][j]-data[1][j])[goodguys]))
          logging.info('data     : {}'.format(np.absolute(data[1][j])[goodguys]))
#      for ok,valA,valB in zip(allclose,data[0],data[1]):
#        logging.debug('{}:\n{}\n{}'.format(ok,valA,valB))

    allclose = True                                                            # start optimistic
    for i in range(1,len(data)):
      allclose &= np.allclose(data[i-1],data[i],rtol,atol)                     # accumulate "pessimism"

    return allclose


  def compare_TableRefCur(self,headingsRef,ref,headingsCur='',cur='',
                               normHeadings='',normType=None,
                               absoluteTolerance=False,perLine=False,skipLines=[]):

    return self.compare_Table(headingsRef,
                              self.fileInReference(ref),
                              headingsRef if headingsCur == '' else headingsCur,
                              self.fileInCurrent(ref if cur == '' else cur),
                              normHeadings,normType,
                              absoluteTolerance,perLine,skipLines)


  def compare_TableCurCur(self,headingsCur0,Cur0,Cur1,
                               headingsCur1='',
                               normHeadings='',normType=None,
                               absoluteTolerance=False,perLine=False,skipLines=[]):

    return self.compare_Table(headingsCur0,
                              self.fileInCurrent(Cur0),
                              headingsCur0 if headingsCur1 == '' else headingsCur1,
                              self.fileInCurrent(Cur1),
                              normHeadings,normType,absoluteTolerance,perLine,skipLines)


  def report_Success(self,culprit):

    ret = culprit

    if culprit == 0:
      count = len(self.variants) if self.options.select is None else len(self.options.select)
      msg = 'Test passed.' if count == 1 else 'All {} tests passed.'.format(count)
    elif culprit == -1:
      msg = 'Warning: could not start test...'
      ret = 0
    else:
      msg = 'Test "{}" failed.'.format(self.variantName(culprit-1))

    logging.critical('\n'.join(['*'*40,msg,'*'*40]) + '\n')
    return ret
