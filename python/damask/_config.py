from io import StringIO
import abc

import yaml

class NiceDumper(yaml.SafeDumper):
    """Make YAML readable for humans."""

    def write_line_break(self, data=None):
        super().write_line_break(data)

        if len(self.indents) == 1:
            super().write_line_break()

    def increase_indent(self, flow=False, indentless=False):
        return super().increase_indent(flow, False)


class Config(dict):
    """YAML-based configuration."""

    def __repr__(self):
        """Show as in file."""
        output = StringIO()
        self.save(output)
        output.seek(0)
        return ''.join(output.readlines())

    @classmethod
    def load(cls,fname):
        """
        Load from yaml file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file for writing.

        """
        try:
            fhandle = open(fname)
        except TypeError:
            fhandle = fname
        return cls(yaml.safe_load(fhandle))

    def save(self,fname,**kwargs):
        """
        Save to yaml file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file for writing.
        **kwargs : dict
            Keyword arguments parsed to yaml.dump.

        """
        try:
            fhandle = open(fname,'w')
        except TypeError:
            fhandle = fname

        if 'width' not in kwargs:
            kwargs['width'] = 256
        if 'default_flow_style' not in kwargs:
            kwargs['default_flow_style'] = None
        fhandle.write(yaml.dump(dict(self),Dumper=NiceDumper,**kwargs))


    @property
    @abc.abstractmethod
    def is_complete(self):
        """Check for completeness."""
        pass

    @property
    @abc.abstractmethod
    def is_valid(self):
        """Check for valid file layout."""
        pass
