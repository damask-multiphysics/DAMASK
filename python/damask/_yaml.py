import copy
from io import StringIO
from collections.abc import Iterable
import abc
from typing import Optional, Union, Any, Type, TypeVar

import numpy as np
import yaml
try:
    from yaml import CSafeLoader as SafeLoader
    from yaml import CSafeDumper as SafeDumper
except ImportError:
    from yaml import SafeLoader                                                                     # type: ignore[assignment]
    from yaml import SafeDumper                                                                     # type: ignore[assignment]

from ._typehints import FileHandle
from . import Rotation
from . import util


MyType = TypeVar('MyType', bound='YAML')

class NiceDumper(SafeDumper):
    """Improve YAML readability for humans."""

    def represent_data(self,
                       data: Any):
        """Cast YAML objects and their subclasses to dict."""
        if isinstance(data, dict) and type(data) is not dict:
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

    def write_line_break(self,
                         data: Optional[str] = None):                                               # not for CSafeDumper
        """From https://github.com/yaml/pyyaml/issues/127."""
        super().write_line_break(data)                                                              # type: ignore[misc]

        if len(self.indents) == 1:                                                                  # type: ignore[attr-defined]
            super().write_line_break()                                                              # type: ignore[misc]

    def increase_indent(self,
                        flow: bool = False,
                        indentless: bool = False):                                                  # not for CSafeDumper
        return super().increase_indent(flow, False)                                                 # type: ignore[misc]


class YAML(dict):
    """YAML-based configuration."""

    def __init__(self,
                 config: Optional[Union[str, dict[str, Any]]] = None,
                 **kwargs):
        """
        New YAML-based configuration.

        Parameters
        ----------
        config : dict or str, optional
            YAML. String needs to be valid YAML.
        **kwargs : arbitrary key–value pairs, optional
            Top-level entries of the configuration.

        Notes
        -----
        Values given as key–value pairs take precedence
        over entries with the same key in 'config'.
        """
        if isinstance(config,str):
            kwargs = yaml.load(config, Loader=SafeLoader) | kwargs
        elif isinstance(config,dict):
            kwargs = config | kwargs

        super().__init__(**kwargs)


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
        other : damask.YAML or dict
            Key–value pairs that update self.

        Returns
        -------
        updated : damask.YAML
            Updated configuration.
        """
        duplicate = self.copy()
        duplicate.update(other)
        return duplicate


    def __ior__(self: MyType,
                other) -> MyType:
        """
        Return self|=other.

        Update configuration with contents of other (in-place).

        Parameters
        ----------
        other : damask.YAML or dict
            Key–value pairs that update self.
        """
        self.update(other)
        return self


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
        updated : damask.YAML
            Updated configuration.
        """
        duplicate = self.copy()
        for k in util.to_list(keys):
            del duplicate[k]
        return duplicate


    @classmethod
    def load(cls: Type[MyType],
             fname: FileHandle) -> MyType:
        """
        Load from YAML file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file to read.

        Returns
        -------
        loaded : damask.YAML
            YAML from file.
        """
        with util.open_text(fname) as fhandle:
            return cls(yaml.load(fhandle, Loader=SafeLoader))


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
        for key,default in [('width',256),
                            ('default_flow_style',None),
                            ('sort_keys',False),
                            ('allow_unicode',True),
                            ('Dumper',NiceDumper)]:
            if key not in kwargs:
                kwargs[key] = default

        with util.open_text(fname,'w') as fhandle:
            fhandle.write(yaml.dump(self,**kwargs))


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
