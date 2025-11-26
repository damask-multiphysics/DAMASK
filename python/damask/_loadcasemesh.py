# SPDX-License-Identifier: AGPL-3.0-or-later
from typing import Optional, Union, Any

from ._typehints import FileHandle
from ._yaml import MaskedMatrixDumper
from . import YAML


class LoadcaseMesh(YAML):
    """Load case for mesh solver."""

    def __init__(self,
                 config: Optional[Union[str,dict[str,Any]]] = None,
                 *,
                 loadstep: Optional[list[dict[str,Any]]] = None):
        """
        New mesh solver load case.

        Parameters
        ----------
        config : dict or str, optional
            Mesh solver load case. String needs to be valid YAML.
        loadstep : list of dict, optional
            Load step configuration.
            Defaults to an empty list if 'config' is not given.
        """
        kwargs: dict[str,Union[dict[str,str],list[dict[str,Any]]]] = {}
        if loadstep is not None:
            kwargs['loadstep'] = loadstep
        elif config is None:
            kwargs['loadstep'] = []

        super().__init__(config,**kwargs)


    def save(self,
             fname: FileHandle,
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
