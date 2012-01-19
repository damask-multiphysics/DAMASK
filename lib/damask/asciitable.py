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
    self.headerLen = 0
    self.info = []
    self.labels = []
    self.data = []

  def output_write(self,
                   what):
    if isinstance(what,list):
      for item in what: self.output_write(item)
    else:
      self.__IO__['output'] += [str(what)]
      self.__IO__['buffered'] or self.output_flush()

  def output_flush(self,
                   clear = True):
    self.__IO__['output'] == [] or self.__IO__['out'].write('\n'.join(self.__IO__['output']) + '\n')
    if clear: self.output_clear()

  def output_clear(self):
    self.__IO__['output'] = []

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
      self.headerLen = int(m.group(1)) + 1
    else:
      self.info      = []
      self.labels    = firstline.split()
      self.headerLen = 1
    self.__IO__['validReadSize'] = len(self.labels)
    self.__IO__['dataStart'] = self.__IO__['in'].tell()

  def head_write(self):
    self.output_write (['%i\theader'%(len(self.info)+1),
                        self.info,
                        '\t'.join(self.labels)])

  def labels_append(self,
                   what):
    if isinstance(what,list):
      for item in what: self.labels_append(item)
    else:               self.labels += [str(what)]

  def info_append(self,
                  what):
    if isinstance(what,list):
      for item in what: self.info_append(item)
    else:               self.info += [str(what)]

  def data_rewind(self):
    self.__IO__['in'].seek(self.__IO__['dataStart'])
    
  def data_read(self):
    line = self.__IO__['in'].readline()
    items = line.split()[:self.__IO__['validReadSize']]                             # get next data row
    self.data = {False:   [],
                  True: items}[len(items) == self.__IO__['validReadSize']]          # take if correct number of entries
    return line != ''

  def data_write(self):
    if isinstance(self.data[0],list):
      self.output_write (['\t'.join(map(str,items)) for items in self.data])
    else:
      self.output_write ('\t'.join(map(str,self.data)))

  def data_append(self,
                  what):
    if isinstance(what,list):
      for item in what: self.data_append(item)
    else:               self.data += [str(what)]
