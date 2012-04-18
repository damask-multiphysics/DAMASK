# $Id$

class ASCIItable():
  '''
     There should be a doc string here  :)
  '''
  import sys
  __slots__ = ['__IO__',
               'info',
               'labels',
               'data',
              ]

# ------------------------------------------------------------------
  def __init__(self,
               fileIn = sys.stdin,
               fileOut = sys.stdout,
               buffered = True):
    self.__IO__ = {'in': fileIn,
                   'out':fileOut,
                   'output':[],
                   'buffered':buffered,
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
    except IOError, e:
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
    if m:
      self.info      = [self.__IO__['in'].readline().strip() for i in xrange(1,int(m.group(1)))]
      self.labels    =  self.__IO__['in'].readline().split()
    else:
      self.info      = []
      self.labels    = firstline.split()
    self.__IO__['validReadSize'] = len(self.labels)
    try:
      self.__IO__['dataStart'] = self.__IO__['in'].tell()
    except IOError:
      pass

# ------------------------------------------------------------------
  def head_write(self):
    return self.output_write (['%i\theader'%(len(self.info)+1),
                              self.info,
                              '\t'.join(self.labels)])

# ------------------------------------------------------------------
  def labels_append(self,
                   what):
    if isinstance(what,list):
      for item in what: self.labels_append(item)
    else:               self.labels += [str(what)]

# ------------------------------------------------------------------
  def info_append(self,
                  what):
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
    line = self.__IO__['in'].readline()
    items = line.split()[:self.__IO__['validReadSize']]                             # get next data row
    self.data = {False:   [],
                  True: items}[len(items) == self.__IO__['validReadSize']]          # take if correct number of entries
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
    except ValueError:
      pass

    return idx
    
# ------------------------------------------------------------------
  def data_clear(self):
    self.data = []

# ------------------------------------------------------------------
  def data_asFloat(self):
    return map(self._transliterateToFloat,self.data)

