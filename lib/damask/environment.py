import os,sys,string,re

class Environment():
  __slots__ = ['pathInfo',
              ]

  def __init__(self,rootRelation = '.'):
    self.pathInfo = {\
                     'acml': '/opt/acml4.4.0',
                     'fftw': '.',
                     'msc':  '/msc',
                    }
    self.get_pathInfo(rootRelation)

  def relPath(self,relative = '.'):
    return os.path.join(self.rootDir(),relative)

  def rootDir(self,rootRelation = '.'):      #getting pathinfo
    damask_root = os.getenv('DAMASK_ROOT')
    if damask_root == '' or damask_root == None: damask_root = os.path.join(os.path.dirname(sys.argv[0]),rootRelation)
    return damask_root

  def binDir(self,rootRelation = '.'):      #getting pathinfo
    damask_bin  = os.getenv('DAMASK_BIN')
    if damask_bin  == '' or damask_bin  == None: damask_bin = self.relPath('bin/')
    return damask_bin
    
  def get_pathInfo(self,rootRelation = '.'):      #getting pathinfo
    damask_root = self.rootDir(rootRelation)

    try:                                                  # check for user-defined pathinfo
      file = open(self.relPath('lib/pathinfo'))
      content = file.readlines()
      file.close()
      for line in content:
        self.pathInfo[line.split()[0].lower()] = os.path.normpath(os.path.join(self.relPath('lib/'),line.split()[1]))
    except:
      pass

