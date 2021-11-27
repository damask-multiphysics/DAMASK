import subprocess
import shlex
import re
from pathlib import Path

_marc_version = 2021.2
_marc_root = '/opt/msc'
_damask_root = str(Path(__file__).parents[3])

class Marc:
    """Wrapper to run DAMASK with MSC Marc."""

    def __init__(self,marc_version=_marc_version,marc_root=_marc_root,damask_root=_damask_root):
        """
        Create a Marc solver object.

        Parameters
        ----------
        version : float
            Marc version

        """
        self.marc_version = marc_version
        self.marc_root    = Path(marc_root)
        self.damask_root = Path(damask_root)

    @property
    def library_path(self):

        path_lib = self.marc_root/f'mentat{self.marc_version}/shlib/linux64'
        if not path_lib.is_dir():
            raise FileNotFoundError(f'library path "{path_lib}" not found')

        return path_lib


    @property
    def tools_path(self):

        path_tools = self.marc_root/f'marc{self.marc_version}/tools'
        if not path_tools.is_dir():
            raise FileNotFoundError(f'tools path "{path_tools}" not found')

        return path_tools


    def submit_job(self, model, job,
                   compile      = False,
                   optimization = ''):
        """
        Assemble command line arguments and call Marc executable.

        Parameters
        ----------
        model : str
            Name of model.
        job : str
            Name of job.
        compile : bool, optional
            Compile DAMASK_Marc user subroutine (and save for future use).
            Defaults to False.
        optimization : str, optional
            Optimization level '' (-O0), 'l' (-O1), or 'h' (-O3).
            Defaults to ''.

        """
        usersub = (self.damask_root/'src/DAMASK_Marc').with_suffix('.f90' if compile else '.marc')
        if not usersub.is_file():
            raise FileNotFoundError(f'subroutine ({"source" if compile else "binary"}) "{usersub}" not found')

        # Define options [see Marc Installation and Operation Guide, pp 23]
        script = f'run_damask_{optimization}mp'

        cmd = f'{self.tools_path/script} -jid {model}_{job} -nprocd 1 -autorst 0 -ci n -cr n -dcoup 0 -b no -v no ' \
            + (f'-u {usersub} -save y' if compile else f'-prog {usersub.with_suffix("")}')

        print(cmd)

        ret = subprocess.run(shlex.split(cmd),capture_output=True)

        try:
            v = int(re.search('Exit number ([0-9]+)',ret.stderr.decode()).group(1))
            if 3004 != v:
                print(ret.stderr.decode())
                print(ret.stdout.decode())
                raise RuntimeError(f'Marc simulation failed ({v})')
        except (AttributeError,ValueError):
            print(ret.stderr.decode())
            print(ret.stdout.decode())
            raise RuntimeError('Marc simulation failed (unknown return value)')
