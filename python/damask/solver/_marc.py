import subprocess
import shlex
import string
from pathlib import Path

import damask

class Marc:
    """Wrapper to run DAMASK with MSCMarc."""

    def __init__(self,version=damask.environment.options['MARC_VERSION']):
        """
        Create a Marc solver object.

        Parameters
        ----------
        version : float
            Marc version

        """
        self.solver  = 'Marc'
        self.version = version

    @property
    def library_path(self):

        path_MSC = damask.environment.options['MSC_ROOT']
        path_lib = Path(f'{path_MSC}/mentat{self.version}/shlib/linux64')

        return path_lib if path_lib.is_dir() else None


    @property
    def tools_path(self):

        path_MSC   = damask.environment.options['MSC_ROOT']
        path_tools = Path(f'{path_MSC}/marc{self.version}/tools')

        return path_tools if path_tools.is_dir() else None


#--------------------------
    def submit_job(self,
                   model,
                   job          = 'job1',
                   logfile      = False,
                   compile      = False,
                   optimization = '',
                  ):


        usersub = damask.environment.root_dir/'src/DAMASK_marc'
        usersub = usersub.parent/(usersub.name + ('.f90' if compile else '.marc'))
        if not usersub.is_file():
            raise FileNotFoundError("DAMASK4Marc ({}) '{}' not found".format(('source' if compile else 'binary'),usersub))

        # Define options [see Marc Installation and Operation Guide, pp 23]
        script = f'run_damask_{optimization}mp'

        cmd = str(self.tools_path/Path(script)) + \
              ' -jid ' + model + '_' + job + \
              ' -nprocd 1  -autorst 0 -ci n  -cr n  -dcoup 0 -b no -v no'

        if compile: cmd += ' -u ' + str(usersub) + ' -save y'
        else:       cmd += ' -prog ' + str(usersub.with_suffix(''))

        print('job submission {} compilation: {}'.format(('with' if compile else 'without'),usersub))
        if logfile: log = open(logfile, 'w')
        print(cmd)
        process = subprocess.Popen(shlex.split(cmd),stdout = log,stderr = subprocess.STDOUT)
        log.close()
        process.wait()

#--------------------------
    def exit_number_from_outFile(self,outFile=None):
        exitnumber = -1
        with open(outFile,'r') as fid_out:
            for line in fid_out:
                if (string.find(line,'tress iteration') != -1):
                    print(line)
                elif (string.find(line,'Exit number')   != -1):
                    substr = line[string.find(line,'Exit number'):len(line)]
                    exitnumber = int(substr[12:16])

        return exitnumber
