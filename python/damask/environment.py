import os

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
        for item in ['DAMASK_NUM_THREADS',
                     'MSC_ROOT',
                     'MARC_VERSION',
                     'ABAQUS_VERSION',
                     ]:
            self.options[item] = os.environ[item] if item in os.environ else None
