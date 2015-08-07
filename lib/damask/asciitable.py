# -*- coding: UTF-8 no BOM -*-

# $Id$

import os,sys
import numpy as np

class ASCIItable():
  '''
     There should be a doc string here  :)
  '''

  __slots__ = ['__IO__',
               'info',
               'labeled',
               'data',
              ]

# ------------------------------------------------------------------
  def __init__(self,
               name    = 'STDIN',
               outname = None,
               buffered  = False,                                                                   # flush writes
               labeled   = True,                                                                    # assume table has labels
               readonly  = False,                                                                   # no reading from file
               writeonly = False,                                                                   # no writing to file
              ):
    self.__IO__ = {'output': [],
                   'buffered':  buffered,
                   'labeled':  labeled,                                                             # header contains labels
                   'labels': [],                                                                    # labels according to file info
                   'readBuffer': [],                                                                # buffer to hold non-advancing reads
                   'dataStart': 0,
                  }

    self.__IO__ .update({'in':  sys.stdin,
                         'out': sys.stdout,
                        } if name == 'STDIN' else
                        {'in':  sys.stdin if writeonly else open(name,'r') ,
                         'out': sys.stdout if readonly else open(outname,'w'),
                        }
                       )
    self.info   = []
    self.labels = []
    self.data   = []

# ------------------------------------------------------------------
  def _transliterateToFloat(self,
                            x):
    try:
      return float(x)
    except:
      return 0.0

# ------------------------------------------------------------------
  def croak(self,
            what, newline = True):
            
    sys.stderr.write(('\n'.join(map(str,what)) if not hasattr(what, "strip")
                                                  and hasattr(what, "__getitem__")
                                                  or  hasattr(what, "__iter__") else str(what))
                    +('\n' if newline else '')),

# ------------------------------------------------------------------
  def close(self,
            dismiss = False):
    self.input_close()
    self.output_flush()
    self.output_close(dismiss)

# ------------------------------------------------------------------
  def input_close(self):
    try:
      self.__IO__['in'].close()
    except:
      pass

# ------------------------------------------------------------------
  def output_write(self,
                   what):
    '''
       aggregate a single row (string) or list of (possibly containing further lists of) rows into output
    '''
    if not isinstance(what, (str, unicode)):
      try:
        for item in what: self.output_write(item)
      except:
        self.__IO__['output'] += [str(what)]
    else:
      self.__IO__['output'] += [what]

    return self.__IO__['buffered'] or self.output_flush()

# ------------------------------------------------------------------
  def output_flush(self,
                   clear = True):
    try:
      self.__IO__['output'] == [] or self.__IO__['out'].write('\n'.join(self.__IO__['output']) + '\n')
    except(IOError) as e:
      return False
    if clear: self.output_clear()
    return True

# ------------------------------------------------------------------
  def output_clear(self):
    self.__IO__['output'] = []

# ------------------------------------------------------------------
  def output_close(self,
                   dismiss = False):
    try:
      self.__IO__['out'].close()
    except:
      pass
    if dismiss and os.path.isfile(self.__IO__['out'].name): os.remove(self.__IO__['out'].name)

# ------------------------------------------------------------------
  def head_read(self):
    '''
       get column labels by either reading
       the first row or, if keyword "head[*]" is present,
       the last line of the header
    '''
    import re

    try:
      self.__IO__['in'].seek(0)
    except:
      pass

    firstline = self.__IO__['in'].readline()
    m = re.search('(\d+)\s+head', firstline.lower())                                                # search for "head" keyword
    if self.__IO__['labeled']:                                                                      # table features labels
      if m:                                                                                         # found header info
        self.info   = [self.__IO__['in'].readline().strip() for i in xrange(1,int(m.group(1)))]
        self.labels = self.__IO__['in'].readline().split()                                          # store labels found in last line
      else:                                                                                         # no header info (but labels)
        self.labels = firstline.split()                                                             # store labels from first line

      self.__IO__['labels'] = list(self.labels)                                                     # backup labels (make COPY, not link)

    else:                                                                                           # no labels present in table
      if m:                                                                                         # found header info
        self.info    = [self.__IO__['in'].readline().strip() for i in xrange(0,int(m.group(1)))]    # all header is info ...
                                                                                                    # ... without any labels
      else:                                                                                         # otherwise file starts with data right away
        try:
          self.__IO__['in'].seek(0)                                                                 # try to rewind
        except:
          self.__IO__['readBuffer'] = firstline                                                     # or at least save data in buffer
    try:
      self.__IO__['dataStart'] = self.__IO__['in'].tell()                                           # current file position is at start of data
    except(IOError):
      pass

# ------------------------------------------------------------------
  def head_write(self,
                 header = True):
    '''
       write current header information (info + labels)
    '''
    
    head = ['{}\theader'.format(len(self.info)+self.__IO__['labeled'])] if header else []
    head.append(self.info)
    if self.__IO__['labeled']: head.append('\t'.join(self.labels))
    
    return self.output_write(head)

# ------------------------------------------------------------------
  def head_getGeom(self):
    '''
       interpret geom header
    '''
    identifiers = {
            'grid':    ['a','b','c'],
            'size':    ['x','y','z'],
            'origin':  ['x','y','z'],
              }
    mappings = {
            'grid':            lambda x: int(x),
            'size':            lambda x: float(x),
            'origin':          lambda x: float(x),
            'homogenization':  lambda x: int(x),
            'microstructures': lambda x: int(x),
              }
    info = {
            'grid':            np.zeros(3,'i'),
            'size':            np.zeros(3,'d'),
            'origin':          np.zeros(3,'d'),
            'homogenization':  0,
            'microstructures': 0,
           }
    extra_header = []

    for header in self.info:
      headitems = map(str.lower,header.split())
      if len(headitems) == 0: continue                                                              # skip blank lines
      if headitems[0] in mappings.keys():
        if headitems[0] in identifiers.keys():
          for i in xrange(len(identifiers[headitems[0]])):
            info[headitems[0]][i] = \
              mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])
        else:
          info[headitems[0]] = mappings[headitems[0]](headitems[1])
      else:
        extra_header.append(header)

    return info,extra_header
    
# ------------------------------------------------------------------
  def labels_append(self,
                    what):
    '''
       add item or list to existing set of labels (and switch on labeling)
    '''
    if not isinstance(what, (str, unicode)):
      try:
        for item in what: self.labels_append(item)
      except:
        self.labels += [str(what)]
    else:
      self.labels += [what]

    self.__IO__['labeled'] = True                                                                  # switch on processing (in particular writing) of labels

# ------------------------------------------------------------------
  def labels_clear(self):
    '''
       delete existing labels and switch to no labeling
    '''
    self.labels = []
    self.__IO__['labeled'] = False

# ------------------------------------------------------------------
  def label_index(self,
                  labels):
    '''
       tell index of column label(s).
       return numpy array if asked for list of labels.
       transparently deals with label positions implicitly given as numbers or their headings given as strings.
    '''
    from collections import Iterable

    if isinstance(labels, Iterable) and not isinstance(labels, str):                                # check whether list of labels is requested
      idx = []
      for label in labels:
        if label != None:
          try:
            idx.append(int(label))                                                                  # column given as integer number?
          except ValueError:
            try:
              idx.append(self.labels.index(label))                                                  # locate string in label list
            except ValueError:
              try:
                idx.append(self.labels.index('1_'+label))                                           # locate '1_'+string in label list
              except ValueError:
               idx.append(-1)                                                                       # not found...
    else:
      try:
        idx = int(labels)
      except ValueError:
        try:
          idx = self.labels.index(labels)
        except ValueError:
          try:
            idx = self.labels.index('1_'+labels)                                                    # locate '1_'+string in label list
          except ValueError:
            idx = None if labels == None else -1

    return np.array(idx) if isinstance(idx,list) else idx

# ------------------------------------------------------------------
  def label_dimension(self,
                      labels):
    '''
       tell dimension (length) of column label(s).
       return numpy array if asked for list of labels.
       transparently deals with label positions implicitly given as numbers or their headings given as strings.
    '''

    from collections import Iterable

    if isinstance(labels, Iterable) and not isinstance(labels, str):                                # check whether list of labels is requested
      dim = []
      for label in labels:
        if label != None:
          myDim = -1
          try:                                                                                      # column given as number?
            idx = int(label)
            myDim = 1                                                                               # if found has at least dimension 1
            if self.labels[idx][:2] == '1_':                                                        # column has multidim indicator?
              while idx+myDim < len(self.labels) and self.labels[idx+myDim][:2] == "%i_"%(myDim+1):
                myDim += 1                                                                          # add while found
          except ValueError:                                                                        # column has string label
            if label in self.labels:                                                                # can be directly found?
              myDim = 1                                                                             # scalar by definition
            elif '1_'+label in self.labels:                                                         # look for first entry of possible multidim object
              idx = self.labels.index('1_'+label)                                                   # get starting column
              myDim = 1                                                                             # (at least) one-dimensional
              while idx+myDim < len(self.labels) and self.labels[idx+myDim][:2] == "%i_"%(myDim+1):
                myDim += 1                                                                          # keep adding while going through object

          dim.append(myDim)
    else:
      dim = -1                                                                                      # assume invalid label
      idx = -1
      try:                                                                                          # column given as number?
        idx = int(labels)
        dim = 1                                                                                     # if found has at least dimension 1
        if self.labels[idx][:2] == '1_':                                                            # column has multidim indicator?
          while idx+dim < len(self.labels) and self.labels[idx+dim][:2] == "%i_"%(dim+1):
            dim += 1                                                                                # add as long as found
      except ValueError:                                                                            # column has string label
        if labels in self.labels:                                                                   # can be directly found?
          dim = 1                                                                                   # scalar by definition
        elif '1_'+labels in self.labels:                                                            # look for first entry of possible multidim object
          idx = self.labels.index('1_'+labels)                                                      # get starting column
          dim = 1                                                                                   # is (at least) one-dimensional
          while idx+dim < len(self.labels) and self.labels[idx+dim][:2] == "%i_"%(dim+1):
            dim += 1                                                                                # keep adding while going through object

    return np.array(dim) if isinstance(dim,list) else dim

# ------------------------------------------------------------------
  def label_indexrange(self,
                       labels):
    '''
       tell index range for given label(s).
       return numpy array if asked for list of labels.
       transparently deals with label positions implicitly given as numbers or their headings given as strings.
    '''

    from collections import Iterable

    start = self.label_index(labels)
    dim   = self.label_dimension(labels)
  
    return map(lambda a,b: xrange(a,a+b), zip(start,dim)) if isinstance(labels, Iterable) and not isinstance(labels, str) \
    else   xrange(start,start+dim)
  
# ------------------------------------------------------------------
  def info_append(self,
                  what):
    '''
       add item or list to existing set of infos
    '''
    if not isinstance(what, (str, unicode)):
      try:
        for item in what: self.info_append(item)
      except:
        self.info += [str(what)]
    else:
      self.info += [what]

# ------------------------------------------------------------------
  def info_clear(self):
    '''
       delete any info block
    '''
    self.info = []

# ------------------------------------------------------------------
  def data_rewind(self):
    self.__IO__['in'].seek(self.__IO__['dataStart'])                                                # position file to start of data section
    self.__IO__['readBuffer'] = []                                                                  # delete any non-advancing data reads
    self.labels = list(self.__IO__['labels'])                                                       # restore label info found in header (as COPY, not link)
    self.__IO__['labeled'] = len(self.labels) > 0

# ------------------------------------------------------------------
  def data_skipLines(self,
                     count):
    '''
       wind forward by count number of lines
    '''
    for i in xrange(count):
      alive = self.data_read()

    return alive

# ------------------------------------------------------------------
  def data_read(self,
                advance = True):
    '''
       read next line (possibly buffered) and parse it into data array
    '''
    if len(self.__IO__['readBuffer']) > 0:
      line = self.__IO__['readBuffer'].pop(0)                                                       # take buffered content
    else:
      line = self.__IO__['in'].readline()                                                           # get next data row from file

    if not advance:
      self.__IO__['readBuffer'].append(line)                                                        # keep line just read in buffer

    if self.__IO__['labeled']:                                                                      # if table has labels
      items = line.split()[:len(self.__IO__['labels'])]                                             # use up to label count (from original file info)
      self.data = items if len(items) == len(self.__IO__['labels']) else []                         # take entries if correct number, i.e. not too few compared to label count
    else:
      self.data = line.split()                                                                      # otherwise take all

    return self.data != []

# ------------------------------------------------------------------
  def data_readArray(self,
                     labels = []):
    '''
       read whole data of all (given) labels as numpy array
    '''

    try:
      self.data_rewind()                                                                            # try to wind back to start of data
    except:
      pass                                                                                          # assume/hope we are at data start already...

    if labels == None or labels == []:
      use = None                                                                                    # use all columns (and keep labels intact)
      labels_missing = []
    else:
      indices    = self.label_index(labels)                                                         # check requested labels ...
      dimensions = self.label_dimension(labels)                                                     # ... and remember their dimension
      present  = np.where(indices >= 0)[0]                                                          # positions in request list of labels that are present ...
      missing  = np.where(indices <  0)[0]                                                          # ... and missing in table
      labels_missing = np.array(labels)[missing]                                                    # labels of missing data

      columns = []
      for i,(c,d) in enumerate(zip(indices[present],dimensions[present])):                          # for all valid labels ...
        columns += range(c,c + \
                          (d if str(c) != str(labels[present[i]]) else \
                           1))                                                                      # ... transparently add all components unless column referenced by number or with explicit dimension
      use = np.array(columns)

      self.labels = list(np.array(self.labels)[use])                                                # update labels with valid subset

    self.data = np.loadtxt(self.__IO__['in'], usecols=use,ndmin=2)

    return labels_missing

# ------------------------------------------------------------------
  def data_write(self,
                 delimiter = '\t'):
    '''
       write current data array and report alive output back
    '''
    if len(self.data) == 0: return True

    if isinstance(self.data[0],list):
      return self.output_write([delimiter.join(map(str,items)) for items in self.data])
    else:
      return self.output_write(delimiter.join(map(str,self.data)))

# ------------------------------------------------------------------
  def data_writeArray(self,
                      format = '%g', delimiter = '\t'):
    '''
       write whole numpy array data
    '''
    return np.savetxt(self.__IO__['out'],self.data,fmt = format,delimiter = delimiter)

# ------------------------------------------------------------------
  def data_append(self,
                  what):
    if not isinstance(what, (str, unicode)):
      try:
        for item in what: self.data_append(item)
      except:
        self.data += [str(what)]
    else:
      self.data += [what]

# ------------------------------------------------------------------
  def data_set(self,
               what, where):
    '''
       update data entry in column "where". grows data array if needed.
    '''
    idx = -1
    try:
      idx = self.label_index(where)
      if len(self.data) <= idx:
        self.data_append(['n/a' for i in xrange(idx+1-len(self.data))])         # grow data if too short
      self.data[idx] = str(what)
    except(ValueError):
      pass

    return idx

# ------------------------------------------------------------------
  def data_clear(self):
    self.data = []

# ------------------------------------------------------------------
  def data_asFloat(self):
    return map(self._transliterateToFloat,self.data)



# ------------------------------------------------------------------
  def microstructure_read(self,
                          grid):
    '''
       read microstructure data (from .geom format)
    '''

    N = grid.prod()                                                                 # expected number of microstructure indices in data
    microstructure = np.zeros(N,'i')                                                # initialize as flat array

    i = 0
    while i < N and self.data_read():
      items = self.data
      if len(items) > 2:
        if   items[1].lower() == 'of': items = [int(items[2])]*int(items[0])
        elif items[1].lower() == 'to': items = range(int(items[0]),1+int(items[2]))
        else:                          items = map(int,items)
      else:                            items = map(int,items)

      s = min(len(items), N-i)                                                      # prevent overflow of microstructure array
      microstructure[i:i+s] = items[:s]
      i += s

    return microstructure
