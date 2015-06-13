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
               'labels',
               'data',
              ]

# ------------------------------------------------------------------
  def __init__(self,
               fileIn = sys.stdin,
               fileOut = sys.stdout,
               buffered = False,                            # flush writes
               labels = True):                              # assume table has labels
    self.__IO__ = {'in': fileIn,
                   'out':fileOut,
                   'output':[],
                   'buffered':buffered,
                   'labels':labels,
                   'validReadSize': 0,
                   'readBuffer': [],                        # buffer to hold non-advancing reads
                   'dataStart': 0,
                  }
    self.info = []
    self.labels = []
    self.data = []


# ------------------------------------------------------------------
  def _transliterateToFloat(self,x):
    try:
      return float(x)
    except:
      return 0.0

# ------------------------------------------------------------------
  def close(self,dismiss = False):
    self.input_close()
    self.output_close(dismiss)

# ------------------------------------------------------------------
  def input_close(self):
    self.__IO__['in'].close()

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
  def output_close(self, dismiss = False):
    self.__IO__['out'].close()
    if dismiss and os.path.isfile(self.__IO__['out'].name): os.remove(self.__IO__['out'].name)

# ------------------------------------------------------------------
  def head_read(self):
    '''
       get column labels by either read the first row, or 
       --if keyword "head[*]" is present-- the last line of the header
    '''
    import re
    try:
      self.__IO__['in'].seek(0)
    except:
      pass
    firstline = self.__IO__['in'].readline()
    m = re.search('(\d+)\s+head', firstline.lower())
    if self.__IO__['labels']:                                                                       # table features labels
      if m:                                                                                         # found header info
        self.info      = [self.__IO__['in'].readline().strip() for i in xrange(1,int(m.group(1)))]
        self.labels    =  self.__IO__['in'].readline().split()
      else:                                                                                         # no header info (but labels)
        self.labels    = firstline.split()
  
      self.__IO__['validReadSize'] = len(self.labels)

    else:                                                                                           # no labels present in table
      if m:                                                                                         # found header info
        self.info      = [self.__IO__['in'].readline().strip() for i in xrange(0,int(m.group(1)))]  # all header is info
                                                                                                    # ... without any labels
    try:
      self.__IO__['dataStart'] = self.__IO__['in'].tell()                                           # current file position is at start of data
    except(IOError):
      pass

    if self.__IO__['validReadSize'] == 0:                                                           # in case no valid data length is known
      self.data_read(advance = False)
      self.__IO__['validReadSize'] = len(self.data)                                                 # assume constant data width from first line

# ------------------------------------------------------------------
  def head_write(self):
    '''
       write current header information (info + labels)
    '''
    if self.__IO__['labels']:
      return self.output_write ([
                                 '%i\theader'%(len(self.info)+1),
                                 self.info,
                                 '\t'.join(self.labels),
                                ])
    else:
      return self.output_write ([
                                 '%i\theader'%(len(self.info)),
                                 self.info,
                                ])

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

    self.__IO__['labels'] = True                                        # switch on processing (in particular writing) of labels

# ------------------------------------------------------------------
  def labels_clear(self):
    '''
       delete existing labels and switch to no labeling
    '''
    self.labels = []
    self.__IO__['labels'] = False

# ------------------------------------------------------------------
  def label_index(self,
                  labels):
    '''
       tell index of column label(s).
       return numpy array if asked for list of labels.
       transparently deals with label positions implicitly given as numbers or their headings given as strings.
    '''
    from collections import Iterable
    
    if isinstance(labels, Iterable) and not isinstance(labels, str):      # check whether list of labels is requested
      idx = []
      for label in labels:
        if label != None:
          try:
            idx.append(int(label))                                      # column given as integer number?
          except ValueError:
            try:
              idx.append(self.labels.index(label))                      # locate string in label list
            except ValueError:
              try:
                idx.append(self.labels.index('1_'+label))               # locate '1_'+string in label list
              except ValueError:
               idx.append(-1)                                           # not found...
    else:
      try:
        idx = int(labels)
      except ValueError:
        try:
          idx = self.labels.index(labels)
        except ValueError:
          try:
            idx = self.labels.index('1_'+labels)                        # locate '1_'+string in label list
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
    
    if isinstance(labels, Iterable) and not isinstance(labels, str):      # check whether list of labels is requested
      dim = []
      for label in labels:
        if label != None:
          myDim = -1
          try:                                                            # column given as number?
            idx = int(label)
            myDim = 1                                                     # if found has at least dimension 1
            if self.labels[idx][:2] == '1_':                              # column has multidim indicator?
              while idx+myDim < len(self.labels) and self.labels[idx+myDim][:2] == "%i_"%(myDim+1):
                myDim += 1                                                # add while found
          except ValueError:                                              # column has string label
            if label in self.labels:                                      # can be directly found?
              myDim = 1                                                   # scalar by definition
            elif '1_'+label in self.labels:                               # look for first entry of possible multidim object
              idx = self.labels.index('1_'+label)                         # get starting column
              myDim = 1                                                   # (at least) one-dimensional
              while idx+myDim < len(self.labels) and self.labels[idx+myDim][:2] == "%i_"%(myDim+1):
                myDim += 1                                                # keep adding while going through object

          dim.append(myDim)
    else:
      dim = -1                                                            # assume invalid label
      idx = -1
      try:                                                                # column given as number?
        idx = int(labels)
        dim = 1                                                           # if found has at least dimension 1
        if self.labels[idx][:2] == '1_':                                  # column has multidim indicator?
          while idx+dim < len(self.labels) and self.labels[idx+dim][:2] == "%i_"%(dim+1):
            dim += 1                                                      # add as long as found
      except ValueError:                                                  # column has string label
        if labels in self.labels:                                         # can be directly found?
          dim = 1                                                         # scalar by definition
        elif '1_'+labels in self.labels:                                  # look for first entry of possible multidim object
          idx = self.labels.index('1_'+labels)                            # get starting column
          dim = 1                                                         # is (at least) one-dimensional
          while idx+dim < len(self.labels) and self.labels[idx+dim][:2] == "%i_"%(dim+1):
            dim += 1                                                      # keep adding while going through object

    return np.array(dim) if isinstance(dim,list) else dim

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
    self.__IO__['in'].seek(self.__IO__['dataStart'])                                  # position file to start of data section
    self.__IO__['readBuffer'] = []                                                    # delete any non-advancing data reads
    
# ------------------------------------------------------------------
  def data_skipLines(self,count):
    '''
       wind forward by count number of lines
    '''
    for i in xrange(count):
      alive = self.data_read()

    return alive

# ------------------------------------------------------------------
  def data_read(self,advance = True):
    '''
       read next line (possibly buffered) and parse it into data array
    '''
    if len(self.__IO__['readBuffer']) > 0:
      line = self.__IO__['readBuffer'].pop(0)                                         # take buffered content
    else:
      line = self.__IO__['in'].readline()                                             # get next data row from file
    
    if not advance:
      self.__IO__['readBuffer'].append(line)                                          # keep line just read in buffer

    if self.__IO__['labels']:
      items = line.split()[:self.__IO__['validReadSize']]                             # use up to valid size (label count)
      self.data = items if len(items) == self.__IO__['validReadSize'] else []         # take if correct number of entries
    else:
      self.data = line.split()                                                        # take all
      
    return self.data != []
    
# ------------------------------------------------------------------
  def data_readLine(self,line):
    '''
       seek beginning of data and wind forward to selected line
    '''
    self.__IO__['in'].seek(self.__IO__['dataStart'])
    for i in xrange(line-1):
      self.__IO__['in'].readline()
    self.data_read()

# ------------------------------------------------------------------
  def data_readArray(self,
                     labels = []):
    '''
       read whole data of all (given) labels as numpy array
    '''

    if not isinstance(labels,list):
      labels = [labels]
    if labels == [None] or labels == []:
      use = np.arange(self.__IO__['validReadSize'])                                                 # use all columns (and keep labels intact)
      labels_missing = []
    else:
      indices = self.label_index(labels)                                                            # check requested labels
      present  = np.where(indices >= 0)[0]                                                          # positions in request list of labels that are present ...
      missing  = np.where(indices <  0)[0]                                                          # ... and missing in table
      labels_missing = np.array(labels)[missing]                                                    # labels of missing data

      columns = []
      for c in indices[present]:                                                                    # for all valid labels ...
        columns += range(c,c+self.label_dimension(c))                                               # ... transparently add all components
      use = np.array(columns)

      self.labels = list(np.array(self.labels)[use])                                                # ... for missing and present columns
      self.__IO__['validReadSize'] = len(use)                                                       # update data width
    
    try:
      self.data_rewind()                                                                            # try to wind back to start of data
    except:
      pass                                                                                          # assume/hope we are at data start already...
    self.data = np.loadtxt(self.__IO__['in'], usecols=use,ndmin=2)
    return labels_missing
    
# ------------------------------------------------------------------
  def data_write(self,delimiter = '\t'):
    '''
       write current data array and report alive output back
    '''
    if len(self.data) == 0: return True
    
    if isinstance(self.data[0],list):
      return self.output_write([delimiter.join(map(str,items)) for items in self.data])
    else:
      return self.output_write(delimiter.join(map(str,self.data)))

# ------------------------------------------------------------------
  def data_writeArray(self,format = '%g',delimiter = '\t'):
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
               what,where):
    idx = -1
    try:
      idx = self.labels.index(where)
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
