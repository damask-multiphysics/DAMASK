# $Id$

class ASCIItable():
  '''
     There should be a doc string here  :)
  '''
  import sys,numpy
  __slots__ = ['__IO__',
               'info',
               'labels',
               'data',
              ]

# ------------------------------------------------------------------
  def __init__(self,
               fileIn = sys.stdin,
               fileOut = sys.stdout,
               buffered = True,
               labels = True):
    self.__IO__ = {'in': fileIn,
                   'out':fileOut,
                   'output':[],
                   'buffered':buffered,
                   'labels':labels,
                   'validReadSize': 0,
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
  def output_write(self,
                   what):
    if isinstance(what,list):
      for item in what: self.output_write(item)
    else:
      self.__IO__['output'] += [str(what)]
      return self.__IO__['buffered'] or self.output_flush()

# ------------------------------------------------------------------
  def output_flush(self,
                   clear = True):
    import sys
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
      self.__IO__['validReadSize'] = len(self.__IO__['in'].readline().split())                      # assume constant data width from first line

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
       add item or list to existing set of labels
    '''
    if isinstance(what,list):
      for item in what: self.labels_append(item)
    else:               self.labels += [str(what)]

    self.__IO__['labels'] = True                                        # switch on processing (in particular writing) of labels

# ------------------------------------------------------------------
  def labels_clear(self):
    self.labels = []

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
          idx.append(self.labels.index(label))
        except(ValueError):
          idx.append(-1)
    else:
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
    if isinstance(what,list):
      for item in what: self.info_append(item)
    else:               self.info += [str(what)]

# ------------------------------------------------------------------
  def info_clear(self):
    self.info = []

# ------------------------------------------------------------------
  def data_rewind(self):
    self.__IO__['in'].seek(self.__IO__['dataStart'])
    
# ------------------------------------------------------------------
  def data_skipLines(self,lines):
    for i in range(lines):
      self.__IO__['in'].readline()

# ------------------------------------------------------------------
  def data_read(self):
    line = self.__IO__['in'].readline()                                               # get next data row
    if self.__IO__['labels']:
      items = line.split()[:self.__IO__['validReadSize']]                             # use up to valid size (label count)
      self.data = {False:   [],
                    True: items}[len(items) == self.__IO__['validReadSize']]          # take if correct number of entries
    else:
      self.data = line.split()                                                        # take all
      
    return self.data != []
    
# ------------------------------------------------------------------
  def data_readLine(self,line):
    self.__IO__['in'].seek(self.__IO__['dataStart'])
    for i in range(line-1):
      self.__IO__['in'].readline()
    self.data_read()
    
# ------------------------------------------------------------------
  def data_write(self):
    if len(self.data) == 0: return
    if isinstance(self.data[0],list):
      return self.output_write (['\t'.join(map(str,items)) for items in self.data])
    else:
      return self.output_write ('\t'.join(map(str,self.data)))

# ------------------------------------------------------------------
  def data_writeArray(self,format):
    import numpy
    '''
       write whole numpy array data
    '''
    return numpy.savetxt(self.__IO__['out'], self.data, fmt=format)

# ------------------------------------------------------------------
  def data_append(self,
                  what):
    if isinstance(what,list):
      for item in what: self.data_append(item)
    else:               self.data += [str(what)]

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

# ------------------------------------------------------------------
  def data_asArray(self,
                   labels = []):
    import numpy
    '''
       read whole data of all (given) labels as numpy array
    '''

    if labels == []: indices = range(self.__IO__['validReadSize'])              # use all columns
    else: indices = self.labels_index(labels)                                   # use specified columns

    self.data_rewind()
    return numpy.loadtxt(self.__IO__['in'], usecols=indices)

