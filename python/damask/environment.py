import os
import re

class Environment():
    __slots__ = [ \
                  'options',
                ]

    def __init__(self):
        """Read and provide values of DAMASK configuration."""
        self.options = {}
        self.get_options()

    def relPath(self,relative = '.'):
        return os.path.join(self.rootDir(),relative)

    def rootDir(self):
        return os.path.normpath(os.path.join(os.path.realpath(__file__),'../../../'))

    def get_options(self):
        with open(self.relPath(self.rootDir()+'/CONFIG')) as configFile:
            for line in configFile:
                line = line.strip()
                line = line[4:] if line.startswith('set ') else line                                   # remove "set" (tcsh) when setting variables
                if line and not line.startswith('#'):
                    items = re.split(r'\s*=\s*',line)
                    if len(items) == 2: 
                        self.options[items[0].upper()] = \
                          re.sub(r'${*DAMASK_ROOT}*',self.rootDir(),os.path.expandvars(items[1]))      # expand all shell variables and DAMASK_ROOT
