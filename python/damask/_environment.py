import os
from pathlib import Path

class Environment:

    def __init__(self):
        """Read and provide values of DAMASK configuration."""
        self.screen_width  = 1024
        self.screen_height =  768
        try:
            import wx
            app = wx.App(False)
            self.screenwidth, self.screenheight = wx.GetDisplaySize()
        except ImportError:
            try:
                import tkinter
                tk = tkinter.Tk()
                self.screen_width  = tk.winfo_screenwidth()
                self.screen_height = tk.winfo_screenheight()
                tk.destroy()
            except:
                pass

    @property
    def options(self):
        options = {}
        for item in ['DAMASK_NUM_THREADS',
                     'MSC_ROOT',
                     'MARC_VERSION',
                     ]:
            options[item] = os.environ[item] if item in os.environ else None

        return options

    @property
    def root_dir(self):
        """Return DAMASK root path."""
        return Path(__file__).parents[2]


    # for compatibility
    def rootDir(self):
        return str(self.root_dir)
