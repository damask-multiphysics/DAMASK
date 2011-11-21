import os,string,re

class DAMASK_TOOLS():

  def check_env(self):
    import os
    if os.getenv('DAMASK_ROOT') is None:
      print('No DAMASK_ROOT environment variable, did you run DAMASK/installation/setup_shellrc?')
      sys.exit(1)
    else:
      return True       

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
      

  def add_data(self, part=None, section=None, data={}):
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
      old_len=len(self.data['crystallite'])
      self.data['crystallite'][label]={'(output)':[[o] for o in output],'__order__':'(output)'}
      if len(self.data['crystallite'])>old_len: # added new label
        self.data['crystallite']['__order__'].append(label)
  
  def add_texture(self, label='',type='', eulers=[], scatter=0., fraction=1.):
      ''' Experimental! Needs expansion to multi-component textures...''' 
      old_len=len(self.data['texture'])
      if type == '(gauss)':
          gauss={type:[['phi1',eulers[0],'Phi',eulers[1], 'phi2',eulers[2],'scatter',scatter,          'fraction',fraction]],'__order__':label}
          self.data['texture'][label]=gauss
          if len(self.data['texture'])>old_len: # added new label
            self.data['texture']['__order__'].append(label)

  def add_phase(self, file='', label='', phase=None):
      ''' USAGE:
          -read phase "label" from file
          OR
          -phase is dict with one key
      '''
      print(file,label,phase)
      old_len=len(self.data['phase'])
      if file and label and (phase is None):
        other=MATERIAL_CONFIG()
        other.read(file=file)
        phase={label:other.data['phase'][label]}
        label=None
        print phase
      if len(phase)==1 and label is None:
        print('Adding phase %s'%phase.keys()[0])
        label=phase.keys()[0]
        self.data['phase'][label]=phase[label]
        if len(self.data['phase'])>old_len: # added new label
          self.data['phase']['__order__'].append(label)
      else: 
        raise Exception('Wrong arguments')

  def add_microstructure(self, label=None,
                               crystallite=None, 
                               phases=None, 
                               textures=None, 
                               fractions=None):
    ''' Experimental! Needs expansion to multi-constituent microstructures...''' 
    old_len=len(self.data['microstructure'])
    c=self.data['crystallite']['__order__'].index(crystallite)+1
    constituent=phases[:]
    for i in range(len(phases)):
      p=self.data['phase']['__order__'].index(phases[i])+1
      t=self.data['texture']['__order__'].index(textures[i])+1           
      f=fractions[i]   
      constituent[i]=['phase','%i'%p,'texture','%i'%t,'fraction','%f'%f]
    self.data['microstructure'][label]={'crystallite':['%i'%c],
               '(constituent)':constituent,
               '__order__':['crystallite','(constituent)']}
    if len(self.data['microstructure'])>old_len: # added new label           
      self.data['microstructure']['__order__'].append(label)
      





    
  