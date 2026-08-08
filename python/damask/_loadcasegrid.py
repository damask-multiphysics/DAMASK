# SPDX-License-Identifier: AGPL-3.0-or-later
from typing import Any

from . import YAML
from ._typehints import FileHandleText
from ._yaml import MaskedMatrixDumper


class LoadcaseGrid(YAML):
    """Load case for grid solver."""

    def __init__(self,
                 config: str | dict[str, Any] | None = None,
                 *,
                 solver: dict[str, str] | None = None,
                 loadstep: list[dict[str, Any]] | None = None):
        """
        New grid solver load case.

        Parameters
        ----------
        config : dict or str, optional
            Grid solver load case. String needs to be valid YAML.
        solver : dict, optional
            Solver configuration.
            Defaults to an empty dict if 'config' is not given.
        loadstep : list of dict, optional
            Load step configuration.
            Defaults to an empty list if 'config' is not given.
        """
        kwargs: dict[str,dict[str,str] | list[dict[str,Any]]] = {}
        default: list | dict
        for arg,value,default in [('solver',solver,{}),('loadstep',loadstep,[])]:                   # type: ignore[assignment]
            if value is not None:
                kwargs[arg] = value
            elif config is None:
                kwargs[arg] = default

        super().__init__(config,**kwargs)


    def save(self,
             fname: FileHandleText,
             **kwargs):
        """
        Save to YAML file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file to write.
        **kwargs : dict
            Keyword arguments parsed to yaml.dump.
        """
        if 'Dumper' not in kwargs:
            kwargs['Dumper'] = MaskedMatrixDumper

        super().save(fname=fname,**kwargs)
