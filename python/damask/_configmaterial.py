import copy

import numpy as np

from . import grid_filters
from . import Config
from . import Lattice
from . import Rotation
from . import Table

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


    @staticmethod
    def load_table(fname,coordinates=None,constituents={},**kwargs):
        """
        Load from an ASCII table.

        Parameters
        ----------
        fname : str, file handle, or damask.Table
            Table that contains material information.
        coordinates : str, optional
            Label of spatial coordiates. Used for sorting and performing a
            sanity check. Default to None, in which case no sorting or checking is
            peformed.
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
        >>> damask.Table.load('small.txt')
            pos  pos  pos   qu   qu    qu    qu   phase    homog
        0    0    0    0  0.19  0.8   0.24 -0.51  Aluminum SX
        1    1    0    0  0.8   0.19  0.24 -0.51  Steel    SX
        >>> cm.load_table('small.txt','pos',{'O':'qu','phase':'phase'},homogenization='h')
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
        t = Table.load(fname)
        if coordinates is not None:
            t = t.sort_by([f'{i}_{coordinates}' for i in range(3,0,-1)])
            grid_filters.coord0_check(t.get(coordinates))
        constituents_ = {k:t.get(v) for k,v in constituents.items()}
        kwargs_       = {k:t.get(v) for k,v in kwargs.items()}

        _,idx = np.unique(np.hstack(list({**constituents_,**kwargs_}.values())),return_index=True,axis=0)

        constituents_ = {k:v[idx].squeeze() for k,v in constituents_.items()}
        kwargs_       = {k:v[idx].squeeze() for k,v in kwargs_.items()}

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
                        Lattice(v['lattice'])
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


    def material_add(self,constituents,**kwargs):
        """
        Add material entries.

        Parameters
        ----------
        constituents : dict
            Entries for 'constituents'. The key is the name and the value specifies
            the label of the data column in the table
        **kwargs
            Keyword arguments where the key is the name and  the value specifies
            the label of the data column in the table

        Examples
        --------
        >>> import damask
        >>> m = damask.ConfigMaterial()
        >>> O = damask.Rotation.from_random(3).as_quaternion()
        >>> phase = ['Aluminum','Steel','Aluminum']
        >>> m.material_add(constituents={'phase':phase,'O':O},homogenization='SX')
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
        c = [{'constituents':u} for u in ConfigMaterial._constituents(**constituents)]
        for k,v in kwargs.items():
            if hasattr(v,'__len__') and not isinstance(v,str):
                if len(v) != len(c):
                    raise ValueError('len mismatch 1')
                for i,vv in enumerate(v):
                    c[i][k] = [w.item() for w in vv] if isinstance(vv,np.ndarray) else vv.item()
            else:
                for i in range(len(c)):
                    c[i][k] = v
        dup = copy.deepcopy(self)
        if 'material' not in dup: dup['material'] = []
        dup['material'] +=c

        return dup


    @staticmethod
    def _constituents(N=1,**kwargs):
        """Construct list of constituents."""
        for v in kwargs.values():
            if hasattr(v,'__len__') and not isinstance(v,str): N_material = len(v)

        if N == 1:
            m = [[{'fraction':1.0}] for _ in range(N_material)]
            for k,v in kwargs.items():
                if hasattr(v,'__len__') and not isinstance(v,str):
                    if len(v) != N_material:
                        raise ValueError('len mismatch 2')
                    for i,vv in enumerate(np.array(v)):
                        m[i][0][k] = [w.item() for w in vv] if isinstance(vv,np.ndarray) else vv.item()
                else:
                    for i in range(N_material):
                        m[i][0][k] = v
            return m
        else:
            raise NotImplementedError
