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
    self.data = {\
              'homogenization': {'__order__': []},
              'microstructure': {'__order__': []},
              'crystallite':    {'__order__': []},
              'phase':          {'__order__': []},
              'texture':        {'__order__': []},
           }
           
  def __repr__(self):
    me = []
    for part in [
                 'homogenization',
                 'microstructure',
                 'crystallite',
                 'phase',
                 'texture',
                ]:
      me += ['','#-----------------------------#','<%s>'%part,'#-----------------------------#',]
      for section in self.data[part]['__order__']:
        me += ['','[%s]\t\t#-----------------------'%section,'',]
        for key in self.data[part][section]['__order__']:
          if key.startswith('(') and key.endswith(')'):                       # multiple (key)
            me += ['%s\t%s'%(key,' '.join(values)) for values in self.data[part][section][key]]
          else:                                                               # plain key
            me += ['%s\t%s'%(key,' '.join(self.data[part][section][key]))]
          
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
                  

  def add_homogenization(self, name=None, type=None):
    
    if type is 'isostrain':
      h={'type':type}
      h['ngrains']=ngrains
      return h
      
  def add_crystallite(self, name=None):
      pass
  def add_texture(self, name=None):
      pass
  def add_phase(self, name=None):
      pass
  def add_microstructure(self, name=None, phases=None, textures=None, fractions=None):
    self.data['microstructure'][name] = {}        # kill existing entry (should warn?)
    

  