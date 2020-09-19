from io import StringIO
import copy

import yaml
import numpy as np

from . import Lattice
from . import Rotation

class NiceDumper(yaml.SafeDumper):
    """Make YAML readable for humans."""

    def write_line_break(self, data=None):
        super().write_line_break(data)

        if len(self.indents) == 1:
            super().write_line_break()

    def increase_indent(self, flow=False, indentless=False):
        return super().increase_indent(flow, False)


class Material(dict):
    """Material configuration."""

    def __repr__(self):
        """Show as in file."""
        output = StringIO()
        self.save(output)
        output.seek(0)
        return ''.join(output.readlines())

    @staticmethod
    def load(fname):
        """Load from yaml file."""
        try:
            fhandle = open(fname)
        except TypeError:
            fhandle = fname
        return Material(yaml.safe_load(fhandle))

    def save(self,fname='material.yaml'):
        """
        Save to yaml file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file for reading.

        """
        try:
            fhandle = open(fname,'w')
        except TypeError:
            fhandle = fname
        fhandle.write(yaml.dump(dict(self),width=256,default_flow_style=None,Dumper=NiceDumper))


    @property
    def is_complete(self):
        """Check for completeness."""
        try:
           ok = len(self['microstructure']) > 0
           for m in self['microstructure']:
               ok &= m['homogenization'] in self['homogenization']
               for c in m['constituents']:
                   c['orientation']
                   ok &= c['phase'] in self['phase']
           for p in self['phase'].values():
               ok &= 'lattice' in p
           return ok
        except KeyError:
            return False


    @property
    def is_valid(self):
        """Check for valid file layout."""
        ok = True

        if 'phase' in self:
            for k,v in zip(self['phase'].keys(),self['phase'].values()):
                if 'lattice' in v:
                    try:
                        Lattice(v['lattice'])
                    except KeyError:
                        s = v['lattice']
                        print(f"Invalid lattice: '{s}' in phase '{k}'")
                        ok = False

        if 'microstructure' in self:
            for i,v in enumerate(self['microstructure']):
                if 'constituents' in v:
                    f = 0.0
                    for c in v['constituents']:
                        f+= float(c['fraction'])
                        if 'orientation' in c:
                            try:
                                Rotation.from_quaternion(c['orientation'])
                            except ValueError:
                                o = c['orientation']
                                print(f"Invalid orientation: '{o}' in microstructure '{i}'")
                                ok = False
                    if not np.isclose(f,1.0):
                        print(f"Invalid total fraction '{f}' in microstructure '{i}'")
                        ok = False

        return ok


    def microstructure_rename_phase(self,mapping,ID=None,constituent=None):
        """
        Change phase name in microstructure.

        Parameters
        ----------
        mapping: dictionary
            Mapping from old name to new name
        ID: list of ints, optional
            Limit renaming to selected microstructure IDs.
        constituent: list of ints, optional
            Limit renaming to selected constituents.

        """
        dup = copy.deepcopy(self)
        for i,m in enumerate(dup['microstructure']):
            if ID and i not in ID: continue
            for c in m['constituents']:
                if constituent is not None and c not in constituent: continue
                try:
                    c['phase'] = mapping[c['phase']]
                except KeyError:
                    continue
        return dup


    def microstructure_rename_homogenization(self,mapping,ID=None):
        """
        Change homogenization name in microstructure.

        Parameters
        ----------
        mapping: dictionary
            Mapping from old name to new name
        ID: list of ints, optional
            Limit renaming to selected homogenization IDs.

        """
        dup = copy.deepcopy(self)
        for i,m in enumerate(dup['microstructure']):
            if ID and i not in ID: continue
            try:
                m['homogenization'] = mapping[m['homogenization']]
            except KeyError:
                continue
        return dup
