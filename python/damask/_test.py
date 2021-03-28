import os
import sys
import shutil
import logging
import logging.config
from optparse import OptionParser
from pathlib import Path

import damask

class Test:
  """
  General class for testing.

  Is sub-classed by the individual tests.
  """

  variants = []

  def __init__(self, **kwargs):
    """New test."""
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

    self.parser = OptionParser(description = f'{self.description} (Test class version: {damask.version})',
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
      self.clean()
      self.prepareAll()

    for variant,object in enumerate(self.variants):
      name = self.variantName(variant)
      if self.options.show:
        logging.critical(f'{variant+1}: {name}')
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
            if self.update(variant) != 0: logging.critical(f'update for "{name}" failed.')
          elif not (self.options.accept or self.compare(variant)):                                    # no update, do comparison
            return variant+1                                                                          # return culprit

        except Exception as e:
          logging.critical(f'exception during variant execution: "{e}"')
          return variant+1                                                                            # return culprit
    return 0

  def clean(self):
    """Delete directory tree containing current results."""
    try:
      shutil.rmtree(self.dirCurrent())
    except FileNotFoundError:
      logging.warning(f'removal of directory "{self.dirCurrent()}" not possible...')

    try:
      os.mkdir(self.dirCurrent())
      return True
    except FileExistsError:
      logging.critical(f'creation of directory "{self.dirCurrent()}" failed.')
      return False

  def prepareAll(self):
    """Do all necessary preparations for the whole test."""
    return True

  def prepare(self,variant):
    """Do all necessary preparations for the run of each test variant."""
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


  def fileInRoot(self,dir,file):
    """Path to a file in the root directory of DAMASK."""
    return str(Path(os.environ['DAMASK_ROOT'])/dir/file)


  def fileInReference(self,file):
    """Path to a file in the refrence directory for the test."""
    return os.path.join(self.dirReference(),file)


  def fileInCurrent(self,file):
    """Path to a file in the current results directory for the test."""
    return os.path.join(self.dirCurrent(),file)


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
      except FileNotFoundError:
        logging.critical(f'error copying {source} to {target}')
        raise FileNotFoundError


  def copy_Reference2Current(self,sourcefiles=[],targetfiles=[]):

    if len(targetfiles) == 0: targetfiles = sourcefiles
    for i,f in enumerate(sourcefiles):
      try:
        shutil.copy2(self.fileInReference(f),self.fileInCurrent(targetfiles[i]))
      except FileNotFoundError:
        logging.critical(f'Reference2Current: Unable to copy file "{f}"')
        raise FileNotFoundError


  def copy_Base2Current(self,sourceDir,sourcefiles=[],targetfiles=[]):

    source = os.path.normpath(os.path.join(self.dirBase,'../../..',sourceDir))
    if len(targetfiles) == 0: targetfiles = sourcefiles
    for i,f in enumerate(sourcefiles):
      try:
        shutil.copy2(os.path.join(source,f),self.fileInCurrent(targetfiles[i]))
      except FileNotFoundError:
        logging.error(os.path.join(source,f))
        logging.critical(f'Base2Current: Unable to copy file "{f}"')
        raise FileNotFoundError


  def copy_Current2Reference(self,sourcefiles=[],targetfiles=[]):

    if len(targetfiles) == 0: targetfiles = sourcefiles
    for i,f in enumerate(sourcefiles):
      try:
        shutil.copy2(self.fileInCurrent(f),self.fileInReference(targetfiles[i]))
      except FileNotFoundError:
        logging.critical(f'Current2Reference: Unable to copy file "{f}"')
        raise FileNotFoundError


  def copy_Current2Current(self,sourcefiles=[],targetfiles=[]):

    for i,f in enumerate(sourcefiles):
      try:
        shutil.copy2(self.fileInReference(f),self.fileInCurrent(targetfiles[i]))
      except FileNotFoundError:
        logging.critical(f'Current2Current: Unable to copy file "{f}"')
        raise FileNotFoundError


  def execute_inCurrentDir(self,cmd,env=None):

    logging.info(cmd)
    out,error = damask.util.execute(cmd,self.dirCurrent())

    logging.info(error)
    logging.debug(out)

    return out,error


  def report_Success(self,culprit):

    ret = culprit

    if culprit == 0:
      count = len(self.variants) if self.options.select is None else len(self.options.select)
      msg = ('Test passed.' if count == 1 else f'All {count} tests passed.') + '\a\a\a'
    elif culprit == -1:
      msg = 'Warning: could not start test...'
      ret = 0
    else:
      msg = f'Test "{self.variantName(culprit-1)}" failed.'

    logging.critical('\n'.join(['*'*40,msg,'*'*40]) + '\n')
    return ret
