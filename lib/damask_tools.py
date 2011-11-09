class DAMASK_TOOLS():
    import os,string
    def check_env(self):
        import os
        if os.getenv('DAMASK_ROOT') is None:
          print('No DAMASK_ROOT environment variable, did you run DAMASK/installation/setup_shellrc?')
          sys.exit(1)
        else:
          return True         