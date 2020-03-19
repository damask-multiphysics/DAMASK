import os
try:
    import tkinter
except:
    tkinter = None

class Environment:

    def __init__(self):
        """Read and provide values of DAMASK configuration."""
        self.options = self._get_options()
        if tkinter:
            tk = tkinter.Tk()
            self.screen_width  = tk.winfo_screenwidth()
            self.screen_height = tk.winfo_screenheight()
            tk.destroy()
        else:
            self.screen_width  = 1024
            self.screen_height =  768

    def relPath(self,relative = '.'):
        """Return absolute path from path relative to DAMASK root."""
        return os.path.join(self.rootDir(),relative)


    def rootDir(self):
        """Return DAMASK root path."""
        return os.path.normpath(os.path.join(os.path.realpath(__file__),'../../../'))


    def _get_options(self):
        options = {}
        for item in ['DAMASK_NUM_THREADS',
                     'MSC_ROOT',
                     'MARC_VERSION',
                     ]:
            options[item] = os.environ[item] if item in os.environ else None

        return options
