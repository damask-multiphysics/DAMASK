import sys
import time
import os
import subprocess
import shlex
from optparse import Option
from queue import Queue
from threading import Thread


import numpy as np

class bcolors:
    """
    ASCII Colors (Blender code).
    
    https://svn.blender.org/svnroot/bf-blender/trunk/blender/build_files/scons/tools/bcolors.py
    http://stackoverflow.com/questions/287871/print-in-terminal-with-colors-using-python
    """

    HEADER    = '\033[95m'
    OKBLUE    = '\033[94m'
    OKGREEN   = '\033[92m'
    WARNING   = '\033[93m'
    FAIL      = '\033[91m'
    ENDC      = '\033[0m'
    BOLD      = '\033[1m'
    DIM       = '\033[2m'
    UNDERLINE = '\033[4m'
    CROSSOUT  = '\033[9m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''
        self.BOLD = ''
        self.UNDERLINE = ''
        self.CROSSOUT = ''
    

# -----------------------------
def srepr(arg,glue = '\n'):
  """Joins arguments as individual lines."""
  if (not hasattr(arg, "strip") and
          (hasattr(arg, "__getitem__") or
           hasattr(arg, "__iter__"))):
     return glue.join(str(x) for x in arg)
  return arg if isinstance(arg,str) else repr(arg)

# -----------------------------
def croak(what, newline = True):
  """Writes formated to stderr."""
  if what is not None:
    sys.stderr.write(srepr(what,glue = '\n') + ('\n' if newline else ''))
  sys.stderr.flush()

# -----------------------------
def report(who = None,
           what = None):
  """Reports script and file name."""
  croak( (emph(who)+': ' if who is not None else '') + (what if what is not None else '') + '\n' )


# -----------------------------
def report_geom(info,
                what = ['grid','size','origin','homogenization','microstructures']):
  """Reports (selected) geometry information."""
  output = {
            'grid'   : 'grid     a b c:  {}'.format(' x '.join(list(map(str,info['grid'  ])))),
            'size'   : 'size     x y z:  {}'.format(' x '.join(list(map(str,info['size'  ])))),
            'origin' : 'origin   x y z:  {}'.format(' : '.join(list(map(str,info['origin'])))),
            'homogenization' :  'homogenization:  {}'.format(info['homogenization']),
            'microstructures' : 'microstructures: {}'.format(info['microstructures']),
           }
  for item in what: croak(output[item.lower()])

# -----------------------------
def emph(what):
  """Formats string with emphasis."""
  return bcolors.BOLD+srepr(what)+bcolors.ENDC

# -----------------------------
def deemph(what):
  """Formats string with deemphasis."""
  return bcolors.DIM+srepr(what)+bcolors.ENDC

# -----------------------------
def delete(what):
  """Formats string as deleted."""
  return bcolors.DIM+srepr(what)+bcolors.ENDC

# -----------------------------
def strikeout(what):
  """Formats string as strikeout."""
  return bcolors.CROSSOUT+srepr(what)+bcolors.ENDC

# -----------------------------
def execute(cmd,
            streamIn = None,
            wd = './',
            env = None):
  """Executes a command in given directory and returns stdout and stderr for optional stdin."""
  initialPath = os.getcwd()
  os.chdir(wd)
  myEnv = os.environ if env is None else env
  process = subprocess.Popen(shlex.split(cmd),
                             stdout = subprocess.PIPE,
                             stderr = subprocess.PIPE,
                             stdin  = subprocess.PIPE,
                             env = myEnv)
  out,error = [i for i in (process.communicate() if streamIn is None
                                                 else process.communicate(streamIn.read().encode('utf-8')))]
  out   = out.decode('utf-8').replace('\x08','')
  error = error.decode('utf-8').replace('\x08','')
  os.chdir(initialPath)
  if process.returncode != 0: raise RuntimeError('{} failed with returncode {}'.format(cmd,process.returncode))
  return out,error

def coordGridAndSize(coordinates):
  """Determines grid count and overall physical size along each dimension of an ordered array of coordinates."""
  dim    = coordinates.shape[1]
  coords = [np.unique(coordinates[:,i]) for i in range(dim)]
  mincorner = np.array(list(map(min,coords)))
  maxcorner = np.array(list(map(max,coords)))
  grid   = np.array(list(map(len,coords)),'i')
  size   = grid/np.maximum(np.ones(dim,'d'), grid-1.0) * (maxcorner-mincorner)                      # size from edge to edge = dim * n/(n-1) 
  size   = np.where(grid > 1, size, min(size[grid > 1]/grid[grid > 1]))                             # spacing for grid==1 equal to smallest among other ones
  delta  = size/grid
  
  N = grid.prod()
  
  if  N != len(coordinates):
    raise ValueError('Data count {} does not match grid {}.'.format(len(coordinates),' x '.join(map(repr,grid))))
    
  if   np.any(np.abs(np.log10((coords[0][1:]-coords[0][:-1])/delta[0])) > 0.01) \
    or np.any(np.abs(np.log10((coords[1][1:]-coords[1][:-1])/delta[1])) > 0.01):
    raise ValueError('regular grid spacing {} violated.'.format(' x '.join(map(repr,delta))))
  if dim==3 and np.any(np.abs(np.log10((coords[2][1:]-coords[2][:-1])/delta[2])) > 0.01):
    raise ValueError('regular grid spacing {} violated.'.format(' x '.join(map(repr,delta))))  
  
  return grid,size
  
# -----------------------------
class extendableOption(Option):
  """
  Used for definition of new option parser action 'extend', which enables to take multiple option arguments.

  Adopted from online tutorial http://docs.python.org/library/optparse.html
  """

  ACTIONS = Option.ACTIONS + ("extend",)
  STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
  TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
  ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

  def take_action(self, action, dest, opt, value, values, parser):
    if action == "extend":
      lvalue = value.split(",")
      values.ensure_value(dest, []).extend(lvalue)
    else:
      Option.take_action(self, action, dest, opt, value, values, parser)

# Print iterations progress
# from https://gist.github.com/aubricus/f91fb55dc6ba5557fbab06119420dd6a
def progressBar(iteration, total, prefix='', bar_length=50):
  """
  Call in a loop to create terminal progress bar.

  @params:
      iteration   - Required  : current iteration (Int)
      total       - Required  : total iterations (Int)
      prefix      - Optional  : prefix string (Str)
      bar_length  - Optional  : character length of bar (Int)
  """
  fraction = iteration / float(total)
  if not hasattr(progressBar, "last_fraction"):                                                     # first call to function
    progressBar.start_time    = time.time()
    progressBar.last_fraction = -1.0
    remaining_time = '   n/a'
  else:
    if fraction <= progressBar.last_fraction or iteration == 0:                                     # reset: called within a new loop
      progressBar.start_time    = time.time()
      progressBar.last_fraction = -1.0
      remaining_time = '   n/a'
    else:
      progressBar.last_fraction = fraction
      remainder = (total - iteration) * (time.time()-progressBar.start_time)/iteration
      remaining_time = '{: 3d}:'.format(int( remainder//3600)) + \
                       '{:02d}:'.format(int((remainder//60)%60)) + \
                       '{:02d}' .format(int( remainder     %60))

  filled_length = int(round(bar_length * fraction))
  bar = '█' * filled_length + '░' * (bar_length - filled_length)

  sys.stderr.write('\r{} {} {}'.format(prefix, bar, remaining_time)),

  if iteration == total: sys.stderr.write('\n')
  sys.stderr.flush()
 

class return_message():
  """Object with formatted return message."""
  
  def __init__(self,message):
    """
    Sets return message.

    Parameters
    ----------
    message : str or list of str
      message for output to screen

    """
    self.message = message
  
  def __repr__(self):
    """Return message suitable for interactive shells."""
    return srepr(self.message)


class ThreadPool:
  """Pool of threads consuming tasks from a queue."""

  class Worker(Thread):
    """Thread executing tasks from a given tasks queue."""

    def __init__(self, tasks):
      """Worker for tasks."""
      Thread.__init__(self)
      self.tasks = tasks
      self.daemon = True
      self.start()

    def run(self):
      while True:
        func, args, kargs = self.tasks.get()
        try:
          func(*args, **kargs)
        except Exception as e:
          # An exception happened in this thread
          print(e)
        finally:
          # Mark this task as done, whether an exception happened or not
          self.tasks.task_done()


  def __init__(self, num_threads):
    """
    Thread pool.

    Parameters
    ----------
    num_threads : int
      number of threads

    """
    self.tasks = Queue(num_threads)
    for _ in range(num_threads):
      self.Worker(self.tasks)

  def add_task(self, func, *args, **kargs):
    """Add a task to the queue."""
    self.tasks.put((func, args, kargs))

  def map(self, func, args_list):
    """Add a list of tasks to the queue."""
    for args in args_list:
      self.add_task(func, args)

  def wait_completion(self):
    """Wait for completion of all the tasks in the queue."""
    self.tasks.join()
