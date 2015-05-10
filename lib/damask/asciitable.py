# -*- coding: UTF-8 no BOM -*-

# $Id$

import sys
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
  def input_close(self):
    return self.__IO__['in'].close()

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
  def output_close(self):
    return self.__IO__['out'].close()

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
    m = re.search('(\d+)\s*head', firstline.lower())
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
  def labels_index(self,
                   labels):
    '''
       tell index of column label(s)
    '''
    if isinstance(labels,list):
      idx = []
      for label in labels:
        try:
          idx.append(label+0)
        except TypeError:
          try:
            idx.append(self.labels.index(label))
          except ValueError:
            idx.append(-1)
    else:
      try:
        idx = labels+0
      except TypeError:
        try:
          idx = self.labels.index(labels)
        except(ValueError):
          idx = -1

    return idx

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

    if labels == []: indices = range(self.__IO__['validReadSize'])                                  # use all columns
    else: 
      indices = self.labels_index(labels)                                                           # use specified columns
      dictionary = dict(zip(indices, labels))
      self.labels_index = range(len(dictionary))
      self.labels = [dictionary[label] for label in sorted(dictionary)]
    try:                   
      self.data_rewind()                                                                            # try to wind back to start of data
    except:
      pass                                                                                          # assume/hope we are at data start already...
    self.data = np.loadtxt(self.__IO__['in'], usecols=indices,ndmin=2)
    return self.data.shape
    
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
