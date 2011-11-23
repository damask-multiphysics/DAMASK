import os,string,re


class DAMASK_TOOLS():

  def check_env(self):
    import os
    if os.getenv('DAMASK_ROOT') is None:
      print('No DAMASK_ROOT environment variable, did you run DAMASK/installation/setup_shellrc?')
      sys.exit(1)
    else:
      return True       


class ASCII_TABLE():
  import sys
  
  __slots__ = ['__IO__',
               'info',
               'labels',
               'data',
              ]


#.............................................
  def __init__(self,
               fileIn = sys.stdin,
               fileOut = sys.stdout,
               buffered = True):

    self.__IO__ = {'in': fileIn,
                   'out':fileOut,
                   'output':[],
                   'buffered':buffered,
                   'validReadSize': 0,
                  }
    self.info = []
    self.labels = []
    self.data = []


#.............................................
  def output_write(self,
                   what):

    if isinstance(what,list):
      for item in what: self.output_write(item)
    else:
      self.__IO__['output'] += [str(what)]
      self.__IO__['buffered'] or self.output_flush()


#.............................................
  def output_flush(self,
                   clear = True):

    self.__IO__['output'] == [] or self.__IO__['out'].write('\n'.join(self.__IO__['output']) + '\n')
    if clear: self.output_clear()


#.............................................
  def output_clear(self):

    self.__IO__['output'] = []


#.............................................
  def head_read(self):
  #  get column labels by either read the first row, or - if keyword "head[*]" is present - the last line of the header

    try:
      self.__IO__['in'].seek(0)
    except:
      pass

    firstline = self.__IO__['in'].readline()
    m = re.search('(\d+)\s*head', firstline.lower())
    if m:
      self.info     = [self.__IO__['in'].readline().strip() for i in xrange(1,int(m.group(1)))]
      self.labels   =  self.__IO__['in'].readline().split()
    else:
      self.info     = []
      self.labels   = firstline.split()
  
    self.__IO__['validReadSize'] = len(self.labels)


#.............................................
  def head_write(self):

    self.output_write (['%i\theader'%(len(self.info)+1),
                        self.info,
                        '\t'.join(self.labels)])


#.............................................
  def labels_append(self,
                   what):

    if isinstance(what,list):
      for item in what: self.labels_append(item)
    else:               self.labels += [str(what)]


#.............................................
  def info_append(self,
                  what):

    if isinstance(what,list):
      for item in what: self.info_append(item)
    else:               self.info += [str(what)]


#.............................................
  def data_read(self):

    line = self.__IO__['in'].readline()
    items = line.split()[:self.__IO__['validReadSize']]                             # get next data row
    self.data = {False:   [],
                  True: items}[len(items) == self.__IO__['validReadSize']]          # take if correct number of entries

    return line != ''

#.............................................
  def data_write(self):

    if isinstance(self.data[0],list):
      self.output_write (['\t'.join(map(str,items)) for items in self.data])
    else:
      self.output_write ('\t'.join(map(str,self.data)))


#.............................................
  def data_append(self,
                  what):

    if isinstance(what,list):
      for item in what: self.data_append(item)
    else:               self.data += [str(what)]





    
class MATERIAL_CONFIG():

  __slots__ = ['data']
  

  def __init__(self):
    self.parts = [
             'homogenization',
             'microstructure',
             'crystallite',
             'phase',
             'texture',
            ]                                       # ordered (!) list of parts
    self.data = {\
              'homogenization': {'__order__': []},
              'microstructure': {'__order__': []},
              'crystallite':    {'__order__': []},
              'phase':          {'__order__': []},
              'texture':        {'__order__': []},
           }
           
  def __repr__(self):
    me = []
    for part in self.parts:
      print 'doing',part
      me += ['','#-----------------------------#','<%s>'%part,'#-----------------------------#',]
      for section in self.data[part]['__order__']:
        me += ['','[%s] %s'%(section,'-'*max(0,27-len(section))),'',]
        for key in self.data[part][section]['__order__']:
          if key.startswith('(') and key.endswith(')'):                       # multiple (key)
            me += ['%s\t%s'%(key,' '.join(values)) for values in self.data[part][section][key]]
          else:                                                               # plain key
            me += ['%s\t%s'%(key,' '.join(map(str,self.data[part][section][key])))]
          
    return '\n'.join(me)

  def parse_data(self, part=None, sections=[], content=None):

    re_part = re.compile(r'^<(.+)>$')                     # pattern for part
    re_sec  = re.compile(r'^\[(.+)\]$')                   # pattern for section

    name_section = ''
    idx_section = 0
    active = False

    for line in content:
      line = line.split('#')[0].strip()                   # kill comments and extra whitespace
      if line:                                            # content survives...
        match_part = re_part.match(line)
        if match_part:                                    # found <part> separator
          active = (match_part.group(1) == part)          # only active in <part>
          continue
        if active:
          match_sec  = re_sec.match(line)
          if match_sec:                                   # found [section]
            name_section = match_sec.group(1)             # remember name ...
            if '__order__' not in self.data[part]: self.data[part]['__order__'] = []
            self.data[part]['__order__'].append(name_section)  # ... and position
            self.data[part][name_section] = {'__order__':[]}
            continue
          
          if sections == [] or name_section in sections:  # respect subset
            items = line.split()
            if items[0] not in self.data[part][name_section]:         # first encounter of key?
              self.data[part][name_section][items[0]] = []            # create item
              self.data[part][name_section]['__order__'].append(items[0])
            if items[0].startswith('(') and items[0].endswith(')'):   # multiple "(key)"
              self.data[part][name_section][items[0]].append(items[1:])
            else:                                                     # plain key
              self.data[part][name_section][items[0]] = items[1:]
                  
  def read(self,file=None):
     f=open(file,'r')
     c=f.readlines()
     f.close()
     for p in self.parts:
       self.parse_data(part=p, content=c)
       
  def write(self,file='material.config', overwrite=False):
     if overwrite is False:
       if os.path.exists(file):
         i=1
         while os.path.exists(file+'_%i'%i):i+=1
         file+='_%i'%i
     print('Writing material data to file %s'%file)
     f=open(file,'w')
     f.write(str(self))
     f.close()
     return file

  def add_data(self, part=None, section=None, data={}):
    '''Generic data adding/updating'''
    if part not in self.parts: raise Exception('invalid part %s'%part)
    if section not in self.data[part]: self.data[part]['__order__'] += [section]
    self.data[part][section] = data
    
  def add_homogenization(self, label='', type='', Ngrains=None):
      if type.lower() == 'isostrain':
          self.add_data(part='homogenization',
                        section=label,
                        data={'type':[type],
                              'Ngrains':[Ngrains],
                              '__order__':['type','Ngrains']
                             }
                       )
      elif type.lower() == 'rgc':
          raise Exception('Please implement me')
      
  def add_crystallite(self, label='', output=[]):
      self.add_data(part='crystallite', 
                    section=label, 
                    data={'(output)':[[o] for o in output],
                          '__order__':'(output)'})      
      
  def add_texture(self, label='',type='', eulers=[], scatter=0., fraction=1.):
      ''' Experimental! Needs expansion to multi-component textures...''' 
      if type == '(gauss)':
          texture={type:[['phi1','%f'%float(eulers[0]),'Phi','%f'%float(eulers[1]), 'phi2','%f'%float(eulers[2]),'scatter','%f'%float(scatter), 'fraction','%f'%fraction]],'__order__':label}
          #self.data['texture'][label]=texture
          #if len(self.data['texture'])>old_len: # added new label
          # self.data['texture']['__order__'].append(label)
      else:
          raise Exception('Please implement me.')      
      self.add_data(part='texture',section=label, data=texture)      

  def add_phase(self, file='', label='', newlabel=None, phase=None):
      ''' USAGE:
           - read phase "label" from file
          OR
           - phase is dict with one key
      ''' 
      if file and label and (phase is None):
        other=MATERIAL_CONFIG()
        other.read(file=file)
        phase={label:other.data['phase'][label]}
        label=None
        print phase
      if len(phase)==1 and label is None:
        if newlabel:
          label=newlabel
        else:  
          label=phase.keys()[0]
        print('Adding phase %s'%label)
        self.add_data(part='phase', section=label, data=phase)
      else: 
        raise Exception('Wrong arguments')

  def add_microstructure(self, label=None,
                               crystallite=None, # label
                               phases=None,      # list of labels
                               textures=None,    # list of labels
                               fractions=None):  # list of floats
    ''' Experimental! Needs expansion to multi-constituent microstructures...''' 
    c=self.data['crystallite']['__order__'].index(crystallite)+1
    constituent=phases[:]
    if fractions is None:
      fractions=[1./len(phases)]*len(phases)
    for i in range(len(phases)):
      p=self.data['phase']['__order__'].index(phases[i])+1
      t=self.data['texture']['__order__'].index(textures[i])+1           
      f=fractions[i]   
      constituent[i]=['phase','%i'%p,'texture','%i'%t,'fraction','%f'%f]
    data={'crystallite':['%i'%c],
               '(constituent)':constituent,
               '__order__':['crystallite','(constituent)']}               
    self.add_data(part='microstructure',section=label,data=data)           
    
    
  def change_value(self,part=None,section=None,key=None,value=None):
    if type(value) is not type([]): 
      if type(value) is not type('s'):
        value='%s'%value
      value=[value]
    newlen=len(value)  
    oldval=self.data[part][section][key]
    oldlen=len(oldval)
    print('changing %s:%s:%s:%s'%(part,section,key,oldval))
    self.data[part][section][key]=value
    print('new: %s'%self.data[part][section][key])
    if newlen is not oldlen:
      print('Length of value was changed from %i to %i!'%(oldlen,newlen))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    