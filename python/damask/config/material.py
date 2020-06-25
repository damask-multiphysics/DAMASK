import re
import os

from damask import util

class Section():
  def __init__(self,data = {'__order__':[]},part = ''):
    """New material.config section."""
    classes = {
                'homogenization':Homogenization,
                'microstructure':Microstructure,
                'crystallite':Crystallite,
                'phase':Phase,
                'texture':Texture,
              }
    self.parameters = {}
    for key in data:
      self.parameters[key] = data[key] if isinstance(data[key], list) else [data[key]]

    if '__order__' not in self.parameters:
      self.parameters['__order__'] = list(self.parameters.keys())
    if part.lower() in classes:
      self.__class__ = classes[part.lower()]
      self.__init__(data)

  def add_multiKey(self,key,data):
    multiKey = '(%s)'%key
    if multiKey not in self.parameters: self.parameters[multiKey] = []
    if multiKey not in self.parameters['__order__']: self.parameters['__order__'] += [multiKey]
    self.parameters[multiKey] += [[item] for item in data] if isinstance(data, list) else [[data]]

  def data(self):
    return self.parameters


class Homogenization(Section):
  def __init__(self,data = {'__order__':[]}):
    """New material.config <homogenization> section."""
    Section.__init__(self,data)


class Crystallite(Section):
  def __init__(self,data = {'__order__':[]}):
    """New material.config <crystallite> section."""
    Section.__init__(self,data)


class Phase(Section):
  def __init__(self,data = {'__order__':[]}):
    """New material.config <Phase> section."""
    Section.__init__(self,data)


class Microstructure(Section):
  def __init__(self,data = {'__order__':[]}):
    """New material.config <microstructure> section."""
    Section.__init__(self,data)


class Texture(Section):
  def __init__(self,data = {'__order__':[]}):
    """New material.config <texture> section."""
    Section.__init__(self,data)

  def add_component(self,theType,properties):

    scatter  = properties['scatter']  if 'scatter'  in list(map(str.lower,list(properties.keys()))) else 0.0
    fraction = properties['fraction'] if 'fraction' in list(map(str.lower,list(properties.keys()))) else 1.0

    try:
      multiKey = theType.lower()
    except AttributeError:
      pass

    if multiKey == 'gauss':
      self.add_multiKey(multiKey,'phi1 %g\tPhi %g\tphi2 %g\tscatter %g\tfraction %g'%(
                                        properties['eulers'][0],
                                        properties['eulers'][1],
                                        properties['eulers'][2],
                                        scatter,
                                        fraction,
                                        )
                        )


class Material():
  """Read, manipulate, and write material.config files."""

  def __init__(self,verbose=True):
    """Generates ordered list of parts."""
    self.parts = [
            'homogenization',
            'crystallite',
            'phase',
            'texture',
            'microstructure',
            ]
    self.data = {
            'homogenization': {'__order__': []},
            'microstructure': {'__order__': []},
            'crystallite':    {'__order__': []},
            'phase':          {'__order__': []},
            'texture':        {'__order__': []},
            }
    self.verbose = verbose

  def __repr__(self):
    """Returns current data structure in material.config format."""
    me = []
    for part in self.parts:
      if self.verbose: print(f'processing <{part}>')
      me += ['',
             '#'*100,
             f'<{part}>',
             '#'*100,
            ]
      for section in self.data[part]['__order__']:
        me += [f'[{section}] {"#"+"-"*max(0,96-len(section))}']
        for key in self.data[part][section]['__order__']:
          if key.startswith('(') and key.endswith(')'):                       # multiple (key)
            me += [f'{key}\t{" ".join(values)}' for values in self.data[part][section][key]]
          else:                                                               # plain key
            me += [f'{key}\t{util.srepr(self.data[part][section][key]," ")}']
    return '\n'.join(me) + '\n'

  def parse(self, part=None, sections=[], content=None):

    re_part = re.compile(r'^<(.+)>$')                     # pattern for part
    re_sec  = re.compile(r'^\[(.+)\]$')                   # pattern for section

    name_section = ''
    active = False

    for line in content:
      line = line.split('#')[0].strip()                   # kill comments and extra whitespace
      line = line.split('/echo/')[0].strip()              # remove '/echo/' tags
      line = line.lower()                                 # be case insensitive

      if line:                                            # content survives...
        match_part = re_part.match(line)
        if match_part:                                    # found <...> separator
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

          if sections == [] or name_section in sections:  # possibly restrict to subset
            items = line.split()
            if items[0] not in self.data[part][name_section]:         # first encounter of key?
              self.data[part][name_section][items[0]] = []            # create item
              self.data[part][name_section]['__order__'].append(items[0])
            if items[0].startswith('(') and items[0].endswith(')'):   # multiple "(key)"
              self.data[part][name_section][items[0]].append(items[1:])
            else:                                                     # plain key
              self.data[part][name_section][items[0]] = items[1:]



  def read(self,filename=None):
    """Read material.config file."""
    def recursiveRead(filename):
      """Takes care of include statements like '{}'."""
      result = []
      re_include  = re.compile(r'^{(.+)}$')
      with open(filename) as f: lines = f.readlines()
      for line in lines:
        match = re_include.match(line.split()[0]) if line.strip() else False
        result += [line] if not match else \
                  recursiveRead(match.group(1) if match.group(1).startswith('/') else
                                os.path.normpath(os.path.join(os.path.dirname(filename),match.group(1))))
      return result

    c = recursiveRead(filename)
    for p in self.parts:
      self.parse(part=p, content=c)

  def write(self,filename='material.config', overwrite=False):
    """Write to material.config."""
    i = 0
    outname = filename
    while os.path.exists(outname) and not overwrite:
     i += 1
     outname = f'{filename}_{i}'

    if self.verbose: print(f'Writing material data to {outname}')
    with open(outname,'w') as f:
      f.write(str(self))
    return outname

  def add_section(self, part=None, section=None, initialData=None, merge=False):
    """Add Update."""
    part    = part.lower()
    section = section.lower()
    if part not in self.parts: raise Exception(f'invalid part {part}')

    if not isinstance(initialData, dict):
      initialData = initialData.data()

    if section not in self.data[part]: self.data[part]['__order__'] += [section]
    if section in self.data[part] and merge:
      for existing in self.data[part][section]['__order__']:                               # replace existing
        if existing in initialData['__order__']:
          if existing.startswith('(') and existing.endswith(')'):                          # multiple (key)
            self.data[part][section][existing] += initialData[existing]                    # add new multiple entries to existing ones
          else:                                                                            # regular key
            self.data[part][section][existing] = initialData[existing]                     # plain replice
      for new in initialData['__order__']:                                                 # merge new content
        if new not in self.data[part][section]['__order__']:
          self.data[part][section][new] = initialData[new]
          self.data[part][section]['__order__'] += [new]
    else:
      self.data[part][section] = initialData




  def add_microstructure(self, section='',
                               components={},      # dict of phase,texture, and fraction lists
                        ):
    """Experimental! Needs expansion to multi-constituent microstructures..."""
    microstructure = Microstructure()
    # make keys lower case (http://stackoverflow.com/questions/764235/dictionary-to-lowercase-in-python)
    components=dict((k.lower(), v) for k,v in components.items())

    for key in ['phase','texture','fraction','crystallite']:
      if isinstance(components[key], list):
        for i, x in enumerate(components[key]):
          try:
            components[key][i] = x.lower()
          except AttributeError:
            pass
      else:
        try:
          components[key] = [components[key].lower()]
        except AttributeError:
          components[key] = [components[key]]

    for (phase,texture,fraction,crystallite) in zip(components['phase'],components['texture'],
                                                    components['fraction'],components['crystallite']):
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
    if not isinstance(value,list):
      if not isinstance(value,str):
        value = '%s'%value
      value = [value]
    newlen = len(value)
    oldval = self.data[part.lower()][section.lower()][key.lower()]
    oldlen = len(oldval)
    print('changing %s:%s:%s from %s to %s '%(part.lower(),section.lower(),key.lower(),oldval,value))
    self.data[part.lower()][section.lower()][key.lower()] = value
    if newlen is not oldlen:
      print('Length of value was changed from %i to %i!'%(oldlen,newlen))


  def add_value(self, part=None,
                  section=None,
                  key=None,
                  value=None):
    if not isinstance(value,list):
      if not isinstance(value,str):
        value = '%s'%value
      value = [value]
    print('adding %s:%s:%s with value %s '%(part.lower(),section.lower(),key.lower(),value))
    self.data[part.lower()][section.lower()][key.lower()] = value
    self.data[part.lower()][section.lower()]['__order__'] += [key.lower()]
