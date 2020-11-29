import os
from pathlib import Path

class Environment:

    @property
    def screen_size(self):
        try:
            import wx
            _ = wx.App(False)                                                                       # noqa
            width, height = wx.GetDisplaySize()
        except ImportError:
            try:
                import tkinter
                tk = tkinter.Tk()
                width  = tk.winfo_screenwidth()
                height = tk.winfo_screenheight()
                tk.destroy()
            except Exception as e:
                width  = 1024
                height =  768

        return (width,height)


    @property
    def options(self):
        options = {}
        for item in ['DAMASK_NUM_THREADS',
                     'MSC_ROOT',
                     'MSC_VERSION',
                     ]:
            options[item] = os.environ[item] if item in os.environ else None

        return options


    @property
    def root_dir(self):
        """Return DAMASK root path."""
        return Path(__file__).parents[2]
