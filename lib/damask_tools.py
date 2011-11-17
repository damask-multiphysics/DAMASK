class DAMASK_TOOLS():
    import os,string
    def check_env(self):
        import os
        if os.getenv('DAMASK_ROOT') is None:
          print('No DAMASK_ROOT environment variable, did you run DAMASK/installation/setup_shellrc?')
          sys.exit(1)
        else:
          return True       


class MATERIAL_CONFIG():
    import os,sys
    
    def __init__(self):
      homogenization=[]
    
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
    def add_microstructure(self, name=None, phase=None, texture=None):
        pass


    
    def write_file(self,label=None):
      fname='material.config'
      if label is not None: fname+='_%s'%label
      #f=open(fname)
      #f.write()
      #f.close()


class HOMOGENIZATION():
    type=None
    