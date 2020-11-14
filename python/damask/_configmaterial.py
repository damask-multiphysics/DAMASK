import copy

import numpy as np

from . import Config
from . import Rotation
from . import Orientation

class ConfigMaterial(Config):
    """Material configuration."""

    def save(self,fname='material.yaml',**kwargs):
        """
        Save to yaml file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path, optional
            Filename or file for writing. Defaults to 'material.yaml'.
        **kwargs
            Keyword arguments parsed to yaml.dump.

        """
        super().save(fname,**kwargs)


    @classmethod
    def load(cls,fname='material.yaml'):
        """
        Load from yaml file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path, optional
            Filename or file for writing. Defaults to 'material.yaml'.

        """
        return super(ConfigMaterial,cls).load(fname)


    @staticmethod
    def from_table(table,constituents={},**kwargs):
        """
        Load from an ASCII table.

        Parameters
        ----------
        table : damask.Table
            Table that contains material information.
        constituents : dict, optional
            Entries for 'constituents'. The key is the name and the value specifies
            the label of the data column in the table
        **kwargs
            Keyword arguments where the key is the name and  the value specifies
            the label of the data column in the table

        Examples
        --------
        >>> import damask
        >>> import damask.ConfigMaterial as cm
        >>> t = damask.Table.load('small.txt')
        >>> t
            pos  pos  pos   qu   qu    qu    qu   phase    homog
        0    0    0    0  0.19  0.8   0.24 -0.51  Aluminum SX
        1    1    0    0  0.8   0.19  0.24 -0.51  Steel    SX
        >>> cm.from_table(t,{'O':'qu','phase':'phase'},homogenization='homog')
        material:
          - constituents:
              - O: [0.19, 0.8, 0.24, -0.51]
                fraction: 1.0
                phase: Aluminum
            homogenization: SX
          - constituents:
              - O: [0.8, 0.19, 0.24, -0.51]
                fraction: 1.0
                phase: Steel
            homogenization: SX

        """
        constituents_ = {k:table.get(v) for k,v in constituents.items()}
        kwargs_       = {k:table.get(v) for k,v in kwargs.items()}

        _,idx = np.unique(np.hstack(list({**constituents_,**kwargs_}.values())),return_index=True,axis=0)

        idx = np.sort(idx)
        constituents_ = {k:np.atleast_1d(v[idx].squeeze()) for k,v in constituents_.items()}
        kwargs_       = {k:np.atleast_1d(v[idx].squeeze()) for k,v in kwargs_.items()}

        return ConfigMaterial().material_add(constituents_,**kwargs_)


    @property
    def is_complete(self):
        """Check for completeness."""
        ok = True
        for top_level in ['homogenization','phase','material']:
            # ToDo: With python 3.8 as prerequisite we can shorten with :=
            ok &= top_level in self
            if top_level not in self: print(f'{top_level} entry missing')

        if ok:
           ok &= len(self['material']) > 0
           if len(self['material']) < 1: print('Incomplete material definition')

        if ok:
            homogenization = set()
            phase          = set()
            for i,v in enumerate(self['material']):
                if 'homogenization' in v:
                    homogenization.add(v['homogenization'])
                else:
                    print(f'No homogenization specified in material {i}')
                    ok = False

                if 'constituents' in v:
                    for ii,vv in enumerate(v['constituents']):
                        if 'O' not in vv:
                            print('No orientation specified in constituent {ii} of material {i}')
                            ok = False
                        if 'phase' in vv:
                            phase.add(vv['phase'])
                        else:
                            print(f'No phase specified in constituent {ii} of material {i}')
                            ok = False

            for k,v in self['phase'].items():
                if 'lattice' not in v:
                    print(f'No lattice specified in phase {k}')
                    ok = False

            for k,v in self['homogenization'].items():
                if 'N_constituents' not in v:
                    print(f'No. of constituents not specified in homogenization {k}')
                    ok = False

            if phase - set(self['phase']):
                print(f'Phase(s) {phase-set(self["phase"])} missing')
                ok = False
            if homogenization - set(self['homogenization']):
                print(f'Homogenization(s) {homogenization-set(self["homogenization"])} missing')
                ok = False

        return ok


    @property
    def is_valid(self):
        """Check for valid file layout."""
        ok = True

        if 'phase' in self:
            for k,v in self['phase'].items():
                if 'lattice' in v:
                    try:
                        Orientation(lattice=v['lattice'])
                    except KeyError:
                        s = v['lattice']
                        print(f"Invalid lattice: '{s}' in phase '{k}'")
                        ok = False

        if 'material' in self:
            for i,v in enumerate(self['material']):
                if 'constituents' in v:
                    f = 0.0
                    for c in v['constituents']:
                        f+= float(c['fraction'])
                        if 'O' in c:
                            try:
                                Rotation.from_quaternion(c['O'])
                            except ValueError:
                                o = c['O']
                                print(f"Invalid orientation: '{o}' in material '{i}'")
                                ok = False
                    if not np.isclose(f,1.0):
                        print(f"Invalid total fraction '{f}' in material '{i}'")
                        ok = False

        return ok


    def material_rename_phase(self,mapping,ID=None,constituent=None):
        """
        Change phase name in material.

        Parameters
        ----------
        mapping: dictionary
            Mapping from old name to new name
        ID: list of ints, optional
            Limit renaming to selected material IDs.
        constituent: list of ints, optional
            Limit renaming to selected constituents.

        """
        dup = copy.deepcopy(self)
        for i,m in enumerate(dup['material']):
            if ID and i not in ID: continue
            for c in m['constituents']:
                if constituent is not None and c not in constituent: continue
                try:
                    c['phase'] = mapping[c['phase']]
                except KeyError:
                    continue
        return dup


    def material_rename_homogenization(self,mapping,ID=None):
        """
        Change homogenization name in material.

        Parameters
        ----------
        mapping: dictionary
            Mapping from old name to new name
        ID: list of ints, optional
            Limit renaming to selected homogenization IDs.

        """
        dup = copy.deepcopy(self)
        for i,m in enumerate(dup['material']):
            if ID and i not in ID: continue
            try:
                m['homogenization'] = mapping[m['homogenization']]
            except KeyError:
                continue
        return dup


    def material_add(self,constituents=None,**kwargs):
        """
        Add material entries.

        Parameters
        ----------
        constituents : dict, optional
            Entries for 'constituents' as key-value pair.
        **kwargs
            Key-value pairs.

        Examples
        --------
        >>> import damask
        >>> O = damask.Rotation.from_random(3).as_quaternion()
        >>> phase = ['Aluminum','Steel','Aluminum']
        >>> m = damask.ConfigMaterial().material_add(constituents={'phase':phase,'O':O},
        ...                                          homogenization='SX')
        >>> m
        material:
          - constituents:
              - O: [0.577764, -0.146299, -0.617669, 0.513010]
                fraction: 1.0
                phase: Aluminum
            homogenization: SX
          - constituents:
              - O: [0.184176, 0.340305, 0.737247, 0.553840]
                fraction: 1.0
                phase: Steel
            homogenization: SX
          - constituents:
              - O: [0.0886257, -0.144848, 0.615674, -0.769487]
                fraction: 1.0
                phase: Aluminum
            homogenization: SX

        """
        length = -1
        for v in kwargs.values():
            if hasattr(v,'__len__') and not isinstance(v,str):
                if length != -1 and len(v) != length:
                    raise ValueError('Cannot add entries of different length')
                else:
                    length = len(v)
        length = max(1,length)

        c = [{} for _ in range(length)] if constituents is None else \
            [{'constituents':u} for u in ConfigMaterial._constituents(**constituents)]
        if len(c) == 1: c = [copy.deepcopy(c[0]) for _ in range(length)]

        if length != 1 and length != len(c):
            raise ValueError('Cannot add entries of different length')

        for k,v in kwargs.items():
            if hasattr(v,'__len__') and not isinstance(v,str):
                for i,vv in enumerate(v):
                    c[i][k] = vv.item() if isinstance(vv,np.generic) else vv
            else:
                for i in range(len(c)):
                    c[i][k] = v
        dup = copy.deepcopy(self)
        dup['material'] = dup['material'] + c if 'material' in dup else c

        return dup


    @staticmethod
    def _constituents(N=1,**kwargs):
        """Construct list of constituents."""
        N_material=1
        for v in kwargs.values():
            if hasattr(v,'__len__') and not isinstance(v,str): N_material = len(v)

        if N == 1:
            m = [[{'fraction':1.0}] for _ in range(N_material)]
            for k,v in kwargs.items():
                if hasattr(v,'__len__') and not isinstance(v,str):
                    if len(v) != N_material:
                        raise ValueError('Cannot add entries of different length')
                    for i,vv in enumerate(np.array(v)):
                        m[i][0][k] = vv.item() if isinstance(vv,np.generic) else vv
                else:
                    for i in range(N_material):
                        m[i][0][k] = v
            return m
        else:
            raise NotImplementedError
