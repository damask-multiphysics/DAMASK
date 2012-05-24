# $Id$
import re

class Section():
  def __init__(self,data = {'__order__':[]},part = ''):
    classes = {
                'homogenization':Homogenization,
                'microstructure':Microstructure,
                'crystallite':Crystallite,
                'phase':Phase,
                'texture':Texture,
              }
    self.parameters = {}
    for key in data:
      if type(data[key]) is not list:
        self.parameters[key] = [data[key]]
      else:
        self.parameters[key] = data[key]

    if '__order__' not in self.parameters:
      self.parameters['__order__'] = self.parameters.keys()
    if part.lower() in classes:
      self.__class__ = classes[part.lower()]
      self.__init__(data)

  def add_multiKey(self,key,data):
    multiKey = '(%s)'%key
    if multiKey not in self.parameters: self.parameters[multiKey] = []
    if multiKey not in self.parameters['__order__']: self.parameters['__order__'] += [multiKey]
    if type(data) == list: self.parameters[multiKey] += [[item] for item in data]
    else:                  self.parameters[multiKey] += [[data]]
      
  def data(self):
    return self.parameters


class Homogenization(Section):
  def __init__(self,data = {'__order__':[]}):
    Section.__init__(self,data)


class Crystallite(Section):
  def __init__(self,data = {'__order__':[]}):
    Section.__init__(self,data)
  

class Phase(Section):
  def __init__(self,data = {'__order__':[]}):
    Section.__init__(self,data)
  

class Microstructure(Section):
  def __init__(self,data = {'__order__':[]}):
    Section.__init__(self,data)
  
  
class Texture(Section):
  def __init__(self,data = {'__order__':[]}):
    Section.__init__(self,data)
    
  def add_component(self,theType,properties):
    
    if 'scatter' not in map(str.lower,properties.keys()):
          scatter = 0.0
    else: scatter = properties['scatter']
    if 'fraction' not in map(str.lower,properties.keys()):
          fraction = 1.0
    else: fraction = properties['fraction']

    multiKey = theType.lower()

    if multiKey == 'gauss':
      self.add_multiKey(multiKey,'phi1 %g\tPhi %g\tphi2 %g\tscatter %g\tfraction %g'%(
                                        properties['eulers'][0],
                                        properties['eulers'][1],
                                        properties['eulers'][2],
                                        scatter,
                                        fraction,
                                        )
                        )

    if multiKey == 'fiber':
      self.add_multiKey(multiKey,'alpha1 %g\talpha2 %g\tbeta1 %g\tbeta2 %g\tscatter %g\tfraction %g'%(
                                        properties['eulers'][0],
                                        properties['eulers'][1],
                                        properties['eulers'][2],
                                        properties['eulers'][3],
                                        scatter,
                                        fraction,
                                        )
                        )

 
 
class Material():
  
  '''
     Reads, manipulates and writes material.config files
  '''
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
    import os
    i = 0
    saveFile = file
    while not overwrite and os.path.exists(saveFile):
     i += 1
     saveFile = file+'_%i'%i

    print('Writing material data to file %s'%saveFile)
    f=open(saveFile,'w')
    f.write(str(self)+'\n')                                          #newline at end
    f.close()
    return saveFile

  def add_section(self, part=None, section=None, object=None, merge = False):
    '''adding/updating'''
    
    if part not in self.parts: raise Exception('invalid part %s'%part)

    if type(object) is dict: data = object
    else: data = object.data()
    
    if section not in self.data[part]: self.data[part]['__order__'] += [section]
    if section in self.data[part] and merge:
      for existing in self.data[part][section]['__order__']:                        # replace existing
        if existing in data['__order__']:
          if existing.startswith('(') and existing.endswith(')'):                   # multiple (key)
            self.data[part][section][existing] += data[existing]                    # add new multiple entries to existing ones
          else:                                                                     # regular key
            self.data[part][section][existing] = data[existing]                     # plain replice
      for new in data['__order__']:                                                 # merge new content
        if new not in self.data[part][section]['__order__']:
          self.data[part][section][new] = data[new]
          self.data[part][section]['__order__'] += [new]
    else:
      self.data[part][section] = data
        
    
      

  def add_microstructure(self, section='',
                               components={},      # dict of phase,texture, and fraction lists
                        ):
    ''' Experimental! Needs expansion to multi-constituent microstructures...''' 
    
    microstructure = Microstructure()
    
    for property in ['phase','texture','fraction','crystallite']:
      if type(components[property]) is not list: components[property] = [components[property]]

    for (phase,texture,fraction,crystallite) in zip(components['phase'],components['texture'],components['fraction'],components['crystallite']):
      microstructure.add_multiKey('constituent','phase %i\ttexture %i\tfraction %g\ncrystallite %i'%(
                                    self.data['phase']['__order__'].index(phase)+1,
                                    self.data['texture']['__order__'].index(texture)+1,
                                    fraction,
                                    self.data['crystallite']['__order__'].index(crystallite)+1))

    self.add_section('microstructure',section,microstructure)

    
  def change_value(self, part=None, 
                         section=None, 
                         key=None, 
                         value=None):
    if type(value) is not type([]): 
      if type(value) is not type('s'):
        value = '%s'%value
      value = [value]
    newlen = len(value)  
    oldval = self.data[part][section][key]
    oldlen = len(oldval)
    print('changing %s:%s:%s:%s'%(part,section,key,oldval))
    self.data[part][section][key] = value
    print('new: %s'%self.data[part][section][key])
    if newlen is not oldlen:
      print('Length of value was changed from %i to %i!'%(oldlen,newlen))
    

    
def ex1():
    mat=Material()
    p=Phase({'constitution':'lump'})   
    t=Texture()    
    t.add_component('gauss',{'eulers':[1,2,3]})
    mat.add_section('phase','phase1',p)
    mat.add_section('texture','tex1',t)
    mat.add_microstructure('mustruct1',{'phase':['phase1']*2,'texture':['tex1']*2,'fraction':[0.2]*2})
    print mat
    mat.write(file='poop')
    mat.write(file='poop',overwrite=True)
    