import copy
from io import StringIO
from collections.abc import Iterable
import abc
import platform
from typing import Optional, Union, Dict, Any, Type, TypeVar

import numpy as np
import yaml
try:
    from yaml import CSafeLoader as SafeLoader
    from yaml import CSafeDumper as SafeDumper
except ImportError:
    from yaml import SafeLoader                                                                     # type: ignore
    from yaml import SafeDumper                                                                     # type: ignore

from ._typehints import FileHandle
from . import Rotation
from . import util

MyType = TypeVar('MyType', bound='Config')

class NiceDumper(SafeDumper):
    """Make YAML readable for humans."""

    def write_line_break(self,
                         data: Optional[str] = None):
        super().write_line_break(data)                                                              # type: ignore

        if len(self.indents) == 1:                                                                  # type: ignore
            super().write_line_break()                                                              # type: ignore

    def increase_indent(self,
                        flow: bool = False,
                        indentless: bool = False):
        return super().increase_indent(flow, False)                                                 # type: ignore

    def represent_data(self,
                       data: Any):
        """Cast Config objects and its subclasses to dict."""
        if isinstance(data, dict) and type(data) != dict:
            return self.represent_data(dict(data))
        if isinstance(data, np.ndarray):
            return self.represent_data(data.tolist())
        if isinstance(data, Rotation):
            return self.represent_data(data.quaternion.tolist())
        if isinstance(data, np.generic):
            return self.represent_data(data.item())

        return super().represent_data(data)

    def ignore_aliases(self,
                       data: Any) -> bool:
        """Do not use references to existing objects."""
        return True

class Config(dict):
    """YAML-based configuration."""

    def __init__(self,
                 config: Optional[Union[str, Dict[str, Any]]] = None,
                 **kwargs):
        """
        New YAML-based configuration.

        Parameters
        ----------
        config : dict or str, optional
            Configuration. String needs to be valid YAML.
        **kwargs: arbitray keyword-value pairs, optional
            Top level entries of the configuration.

        Notes
        -----
        Values given as keyword-value pairs take precedence
        over entries with the same keyword in 'config'.

        """
        if int(platform.python_version_tuple()[1]) >= 9:
            if isinstance(config,str):
                kwargs = yaml.load(config, Loader=SafeLoader) | kwargs
            elif isinstance(config,dict):
                kwargs = config | kwargs                                                            # type: ignore

            super().__init__(**kwargs)
        else:
            if isinstance(config,str):
                c = yaml.load(config, Loader=SafeLoader)
            elif isinstance(config,dict):
                c = config.copy()
            else:
                c = {}
            c.update(kwargs)

            super().__init__(**c)


    def __repr__(self) -> str:
        """
        Return repr(self).

        Show as in file.

        """
        output = StringIO()
        self.save(output)
        output.seek(0)
        return ''.join(output.readlines())


    def __copy__(self: MyType) -> MyType:
        """
        Return deepcopy(self).

        Create deep copy.

        """
        return copy.deepcopy(self)

    copy = __copy__


    def __or__(self: MyType,
               other) -> MyType:
        """
        Return self|other.

        Update configuration with contents of other.

        Parameters
        ----------
        other : damask.Config or dict
            Key-value pairs that update self.

        Returns
        -------
        updated : damask.Config
            Updated configuration.

        Note
        ----
        This functionality is a backport for Python 3.8

        """
        duplicate = self.copy()
        duplicate.update(other)
        return duplicate


    def __ior__(self: MyType,
                other) -> MyType:
        """
        Return self|=other.

        Update configuration with contents of other (in-place).

        """
        return self.__or__(other)


    def delete(self: MyType,
               keys: Union[Iterable, str]) -> MyType:
        """
        Remove configuration keys.

        Parameters
        ----------
        keys : iterable or scalar
            Label of the key(s) to remove.

        Returns
        -------
        updated : damask.Config
            Updated configuration.

        """
        duplicate = self.copy()
        for k in keys if isinstance(keys, Iterable) and not isinstance(keys, str) else [keys]:
            del duplicate[k]
        return duplicate


    @classmethod
    def load(cls: Type[MyType],
             fname: FileHandle) -> MyType:
        """
        Load from yaml file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file for writing.

        Returns
        -------
        loaded : damask.Config
            Configuration from file.

        """
        return cls(yaml.load(util.open_text(fname), Loader=SafeLoader))

    def save(self,
             fname: FileHandle,
             **kwargs):
        """
        Save to yaml file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file for writing.
        **kwargs : dict
            Keyword arguments parsed to yaml.dump.

        """
        if 'width' not in kwargs:
            kwargs['width'] = 256
        if 'default_flow_style' not in kwargs:
            kwargs['default_flow_style'] = None
        if 'sort_keys' not in kwargs:
            kwargs['sort_keys'] = False

        fhandle = util.open_text(fname,'w')
        try:
            fhandle.write(yaml.dump(self,Dumper=NiceDumper,**kwargs))
        except TypeError:                                                                           # compatibility with old pyyaml
            del kwargs['sort_keys']
            fhandle.write(yaml.dump(self,Dumper=NiceDumper,**kwargs))


    @property
    @abc.abstractmethod
    def is_complete(self):
        """Check for completeness."""
        raise NotImplementedError


    @property
    @abc.abstractmethod
    def is_valid(self):
        """Check for valid file layout."""
        raise NotImplementedError
