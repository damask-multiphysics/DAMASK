import os
import sys
import re
import shlex
from collections.abc import Iterable

import numpy as np

# ------------------------------------------------------------------
class ASCIItable():
  """Read and write to ASCII tables."""

  tmpext = '_tmp'                                                                                   # filename extension for in-place access

# ------------------------------------------------------------------
  def __init__(self,
               name,
               labeled  = True,                                                                     # assume table has labels
               readonly = False,                                                                    # no reading from file
              ):
    """Read and write to ASCII tables."""
    self.__IO__ = {'output': [],
                   'labeled':  labeled,                                                             # header contains labels
                   'tags': [],                                                                      # labels according to file info
                   'dataStart': 0,
                  }

    self.__IO__['inPlace'] = name and not readonly
    outname = name + self.tmpext if self.__IO__['inPlace'] else None                                # transparently create tmp file

    try:
      self.__IO__['in'] = (open(   name,'r') if os.access(   name, os.R_OK) else None) if name else sys.stdin
    except TypeError:
      self.__IO__['in'] = name

    try:
      self.__IO__['out'] = (open(outname,'w') if (not os.path.isfile(outname) or
                                                      os.access(     outname, os.W_OK)
                                                 ) and
                                                 (not self.__IO__['inPlace'] or
                                                  not os.path.isfile(name)   or
                                                      os.access(     name, os.W_OK)
                                                 ) else None) if outname else sys.stdout
    except TypeError:
      self.__IO__['out'] = outname

    self.info  = []
    self.tags  = []
    self.data  = []
    self.line  = ''

    if   self.__IO__['in']  is None \
      or self.__IO__['out'] is None: raise IOError                                                 # complain if any required file access not possible


# ------------------------------------------------------------------
  def _removeCRLF(self,
                  string):
    """Delete any carriage return and line feed from string."""
    try:
      return string.replace('\n','').replace('\r','')
    except AttributeError:
      return str(string)

# ------------------------------------------------------------------
  def _quote(self,
             what):
    """Quote empty or white space-containing output."""
    return '{quote}{content}{quote}'.format(
             quote   = ('"' if str(what)=='' or re.search(r"\s",str(what)) else ''),
             content = what)

# ------------------------------------------------------------------
  def close(self,
            dismiss = False):
    if self.__IO__['in'] != sys.stdin: self.__IO__['in'].close()
    self.output_flush()
    if self.__IO__['out'] != sys.stdout: self.__IO__['out'].close()
    if dismiss and os.path.isfile(self.__IO__['out'].name):
      os.remove(self.__IO__['out'].name)
    elif self.__IO__['inPlace']:
      os.rename(self.__IO__['out'].name, self.__IO__['out'].name[:-len(self.tmpext)])

# ------------------------------------------------------------------
  def output_write(self,
                   what):
    """Aggregate a single row (string) or list of (possibly containing further lists of) rows into output."""
    if isinstance(what, str):
      self.__IO__['output'] += [what]
    else:
      try:
        for item in what: self.output_write(item)
      except TypeError:
        self.__IO__['output'] += [str(what)]

    return self.output_flush()

# ------------------------------------------------------------------
  def output_flush(self,
                   clear = True):
    try:
      self.__IO__['output'] == [] or self.__IO__['out'].write('\n'.join(self.__IO__['output']) + '\n')
    except IOError:
      return False
    if clear: self.__IO__['output'] = []
    return True

# ------------------------------------------------------------------
  def head_read(self):
    """
    Get column labels.

    by either reading the first row or,
    if keyword "head[*]" is present, the last line of the header
    """
    try:
      self.__IO__['in'].seek(0)
    except IOError:
      pass

    firstline = self.__IO__['in'].readline().strip()
    m = re.search(r'(\d+)\s+head', firstline.lower())                                                # search for "head" keyword

    if m:                                                                                           # proper ASCIItable format

      if self.__IO__['labeled']:                                                                    # table features labels

        self.info   = [self.__IO__['in'].readline().strip() for i in range(1,int(m.group(1)))]
        self.tags = shlex.split(self.__IO__['in'].readline())                                       # store tags found in last line

      else:

        self.info    = [self.__IO__['in'].readline().strip() for i in range(0,int(m.group(1)))]     # all header is info ...

    else:                                                                                           # other table format
      self.__IO__['in'].seek(0)

      while self.data_read(advance = False, respectLabels = False):
        if self.line[0] in ['#','!','%','/','|','*','$']:                                           # "typical" comment indicators
          self.info_append(self.line)                                                               # store comment as info
          self.data_read()                                                                          # wind forward one line
        else: break                                                                                 # last line of comments

      if self.__IO__['labeled']:                                                                    # table features labels
        self.tags = self.data                                                                       # get tags from last line in "header"...
        self.data_read()                                                                            # ...and remove from buffer

    if self.__IO__['labeled']:                                                                      # table features tags
      self.__IO__['tags'] = list(self.tags)                                                         # backup tags (make COPY, not link)

    try:
      self.__IO__['dataStart'] = self.__IO__['in'].tell()                                           # current file position is at start of data
    except IOError:
      pass

# ------------------------------------------------------------------
  def head_write(self,
                 header = True):
    """Write current header information (info + labels)."""
    head = [f"{len(self.info)+self.__IO__['labeled']}\theader"] if header else []
    head.append(self.info)
    if self.__IO__['labeled']:
      head.append('\t'.join(map(self._quote,self.tags)))
      if len(self.tags) == 0: raise ValueError('no labels present.')

    return self.output_write(head)

# ------------------------------------------------------------------
  def labels_append(self,
                    what,
                    reset = False):
    """Add item or list to existing set of labels (and switch on labeling)."""
    if isinstance(what, str):
      self.tags += [self._removeCRLF(what)]
    else:
      try:
        for item in what: self.labels_append(item)
      except TypeError:
        self.tags += [self._removeCRLF(str(what))]

    self.__IO__['labeled'] = True                                                                  # switch on processing (in particular writing) of tags
    if reset: self.__IO__['tags'] = list(self.tags)                                                # subsequent data_read uses current tags as data size

# ------------------------------------------------------------------
  def labels_clear(self):
    """Delete existing labels and switch to no labeling."""
    self.tags = []
    self.__IO__['labeled'] = False

# ------------------------------------------------------------------
  def labels(self,
             tags = None,
             raw = False):
    """
    Tell abstract labels.

    "x" for "1_x","2_x",... unless raw output is requested.
    operates on object tags or given list.
    """
    if tags is None: tags = self.tags

    if isinstance(tags, Iterable) and not raw:                                                    # check whether list of tags is requested
      id = 0
      dim = 1
      labelList = []

      while id < len(tags):
        if not tags[id].startswith('1_'):
          labelList.append(tags[id])
        else:
          label = tags[id][2:]                                                                    # get label
          while id < len(tags) and tags[id] == f'{dim}_{label}':                                  # check successors
            id  += 1                                                                              # next label...
            dim += 1                                                                              # ...should be one higher dimension
          labelList.append(label)                                                                 # reached end --> store
          id -= 1                                                                                 # rewind one to consider again

        id += 1
        dim = 1

    else:
      labelList = self.tags

    return labelList

# ------------------------------------------------------------------
  def label_index(self,
                  labels):
    """
    Tell index of column label(s).

    return numpy array if asked for list of labels.
    transparently deals with label positions implicitly given as numbers or their headings given as strings.
    """
    if isinstance(labels, Iterable) and not isinstance(labels, str):                                # check whether list of labels is requested
      idx = []
      for label in labels:
        if label is not None:
          try:
            idx.append(int(label)-1)                                                                # column given as integer number?
          except ValueError:
            label = label[1:-1] if label[0] == label[-1] and label[0] in ('"',"'") else label       # remove outermost quotations
            try:
              idx.append(self.tags.index(label))                                                    # locate string in label list
            except ValueError:
              try:
                idx.append(self.tags.index('1_'+label))                                             # locate '1_'+string in label list
              except ValueError:
               idx.append(-1)                                                                       # not found...
    else:
      try:
        idx = int(labels)-1                                                                         # offset for python array indexing
      except ValueError:
        try:
          labels = labels[1:-1] if labels[0] == labels[-1] and labels[0] in ('"',"'") else labels   # remove outermost quotations
          idx = self.tags.index(labels)
        except ValueError:
          try:
            idx = self.tags.index('1_'+labels)                                                      # locate '1_'+string in label list
          except ValueError:
            idx = None if labels is None else -1

    return np.array(idx) if isinstance(idx,Iterable) else idx

# ------------------------------------------------------------------
  def label_dimension(self,
                      labels):
    """
    Tell dimension (length) of column label(s).

    return numpy array if asked for list of labels.
    transparently deals with label positions implicitly given as numbers or their headings given as strings.
    """
    listOfLabels = isinstance(labels, Iterable) and not isinstance(labels, str)                   # check whether list of labels is requested
    if not listOfLabels: labels = [labels]

    dim = []
    for label in labels:
      if label is not None:
        myDim = -1
        try:                                                                                      # column given as number?
          idx = int(label)-1
          myDim = 1                                                                               # if found treat as single column of dimension 1
        except ValueError:                                                                        # column has string label
          label = label[1:-1] if label[0] == label[-1] and label[0] in ('"',"'") else label       # remove outermost quotations
          if label in self.tags:                                                                  # can be directly found?
            myDim = 1                                                                             # scalar by definition
          elif '1_'+label in self.tags:                                                           # look for first entry of possible multidim object
            idx = self.tags.index('1_'+label)                                                     # get starting column
            myDim = 1                                                                             # (at least) one-dimensional
            while idx+myDim < len(self.tags) and self.tags[idx+myDim].startswith("%i_"%(myDim+1)):
              myDim += 1                                                                          # keep adding while going through object

        dim.append(myDim)

    return np.array(dim) if listOfLabels else dim[0]

# ------------------------------------------------------------------
  def label_indexrange(self,
                       labels):
    """
    Tell index range for given label(s).

    return numpy array if asked for list of labels.
    transparently deals with label positions implicitly given as numbers or their headings given as strings.
    """
    start = self.label_index(labels)
    dim   = self.label_dimension(labels)

    return np.hstack([range(s,s+d) for s,d in zip(start,dim)]).astype(int) \
        if isinstance(labels, Iterable) and not isinstance(labels, str)    \
      else range(start,start+dim)

# ------------------------------------------------------------------
  def info_append(self,
                  what):
    """Add item or list to existing set of infos."""
    if isinstance(what, str):
      self.info += [self._removeCRLF(what)]
    else:
      try:
        for item in what: self.info_append(item)
      except TypeError:
        self.info += [self._removeCRLF(str(what))]

# ------------------------------------------------------------------
  def info_clear(self):
    """Delete any info block."""
    self.info = []

# ------------------------------------------------------------------
  def data_rewind(self):
    self.__IO__['in'].seek(self.__IO__['dataStart'])                                                # position file to start of data section
    self.tags = list(self.__IO__['tags'])                                                           # restore label info found in header (as COPY, not link)
    self.__IO__['labeled'] = len(self.tags) > 0

# ------------------------------------------------------------------
  def data_read(self,
                advance = True,
                respectLabels = True):
    """Read next line and parse it into data array."""
    self.line = self.__IO__['in'].readline().strip()

    self.line = self.line.rstrip('\n')

    if self.__IO__['labeled'] and respectLabels:                                                    # if table has labels
      items = shlex.split(self.line)[:len(self.__IO__['tags'])]                                     # use up to label count (from original file info)
      self.data = items if len(items) == len(self.__IO__['tags']) else []                           # take entries if label count matches
    else:
      self.data = shlex.split(self.line)                                                            # otherwise take all

    return self.data != []

# ------------------------------------------------------------------
  def data_readArray(self,
                     labels = []):
    """Read whole data of all (given) labels as numpy array."""
    try:
      self.data_rewind()                                                                            # try to wind back to start of data
    except IOError:
      pass                                                                                          # assume/hope we are at data start already...

    if labels is None or labels == []:
      use = None                                                                                    # use all columns (and keep labels intact)
      labels_missing = []
    else:
      if isinstance(labels, str) or not isinstance(labels, Iterable):                               # check whether labels are a list or single item
        labels = [labels]
      indices    = self.label_index(labels)                                                         # check requested labels ...
      dimensions = self.label_dimension(labels)                                                     # ... and remember their dimension
      present  = np.where(indices >= 0)[0]                                                          # positions in request list of labels that are present ...
      missing  = np.where(indices <  0)[0]                                                          # ... and missing in table
      labels_missing = np.array(labels)[missing]                                                    # labels of missing data

      columns = []
      for i,(c,d) in enumerate(zip(indices[present],dimensions[present])):                          # for all valid labels ...
        # ... transparently add all components unless column referenced by number or with explicit dimension
        columns += list(range(c,c +
                          (d if str(c) != str(labels[present[i]]) else
                           1)))
      use = np.array(columns) if len(columns) > 0 else None

      self.tags = list(np.array(self.__IO__['tags'])[use])                                         # update labels with valid subset

    self.data = np.loadtxt(self.__IO__['in'],usecols=use,ndmin=2)

    return labels_missing

# ------------------------------------------------------------------
  def data_write(self):
    """Write current data array and report alive output back."""
    if len(self.data) == 0: return True

    if isinstance(self.data[0],list):
      return self.output_write(['\t'.join(map(self._quote,items)) for items in self.data])
    else:
      return self.output_write( '\t'.join(map(self._quote,self.data)))

# ------------------------------------------------------------------
  def data_writeArray(self):
    """Write whole numpy array data."""
    for row in self.data:
      try:
        output = list(map(repr,row))
      except Exception:
        output = [repr(row)]

      try:
        self.__IO__['out'].write('\t'.join(output) + '\n')
      except Exception:
        pass

# ------------------------------------------------------------------
  def data_append(self,
                  what):
    if isinstance(what, str):
      self.data += [what]
    else:
      try:
        for item in what: self.data_append(item)
      except TypeError:
        self.data += [str(what)]
