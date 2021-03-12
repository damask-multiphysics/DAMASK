from os import path

import numpy as np
import h5py

from . import Config
from . import Rotation
from . import Orientation
from . import util

class ConfigMaterial(Config):
    """Material configuration."""

    _defaults = {'material': [],
                 'homogenization': {},
                 'phase': {}}

    def __init__(self,d=_defaults):
        """Initialize object with default dictionary keys."""
        super().__init__(d)


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
    def from_table(table,**kwargs):
        """
        Load from an ASCII table.

        Parameters
        ----------
        table : damask.Table
            Table that contains material information.
        **kwargs
            Keyword arguments where the key is the name and the value specifies
            the label of the data column in the table.

        Examples
        --------
        >>> import damask
        >>> import damask.ConfigMaterial as cm
        >>> t = damask.Table.load('small.txt')
        >>> t
            pos  pos  pos   qu   qu    qu    qu   phase    homog
        0    0    0    0  0.19  0.8   0.24 -0.51  Aluminum SX
        1    1    0    0  0.8   0.19  0.24 -0.51  Steel    SX
        1    1    1    0  0.8   0.19  0.24 -0.51  Steel    SX
        >>> cm.from_table(t,O='qu',phase='phase',homogenization='homog')
        material:
          - constituents:
              - O: [0.19, 0.8, 0.24, -0.51]
                v: 1.0
                phase: Aluminum
            homogenization: SX
          - constituents:
              - O: [0.8, 0.19, 0.24, -0.51]
                v: 1.0
                phase: Steel
            homogenization: SX
        homogenization: {}
        phase: {}

        """
        kwargs_       = {k:table.get(v) for k,v in kwargs.items()}

        _,idx = np.unique(np.hstack(list(kwargs_.values())),return_index=True,axis=0)
        idx = np.sort(idx)
        kwargs_ = {k:np.atleast_1d(v[idx].squeeze()) for k,v in kwargs_.items()}

        return ConfigMaterial().material_add(**kwargs_)


    @staticmethod
    def load_DREAM3D(fname,base_group,data_group,ori_data,phase_id,phase_name):
        """
        Load material data from DREAM3D file.

        The parts of homogenization and phase need to be added by the user.

        Parameters
        ----------
        fname : str
            path to the DREAM3D file.
        base_group : str
            Name of the group (folder) below 'DataContainers',
            for example 'SyntheticVolumeDataContainer'.
        data_group : str
            Name of the group (folder) having relevant data for conversion,
            for example 'Grain Data' or 'CellData'.
        ori_data : str
            Name of the dataset having orientation data (working with Euler Angles in dream3D file),
            For example 'EulerAngles'.
        phase_id : str
            Name of the dataset containing phase IDs for each grain,
            for example 'Phases'.
        phase_name : list
            List with name of the phases.

        Examples
        --------
        for grain based data with single phase
        >>> import damask
        >>> import damask.ConfigMaterial as cm
        >>> cm.load_from_Dream3D('20grains16x16x16.dream3D','SyntheticVolumeDataContainer', 'Grain Data',
        ...                      'EulerAngles','Phases',['Ferrite'])

        for point based data with single phase
        >>> import damask
        >>> import damask.ConfigMaterial as cm
        >>> cm.load_from_Dream3D('20grains16x16x16.dream3D','SyntheticVolumeDataContainer', 'CellData',
        ...                       'EulerAngles','Phases',['Ferrite'])

        for grain based data with dual phase
        >>> import damask
        >>> import damask.ConfigMaterial as cm
        >>> cm.load_from_Dream3D('20grains16x16x16.dream3D','SyntheticVolumeDataContainer', 'Grain Data',
        ...                       'EulerAngles','Phases',['Ferrite','Martensite'])

        for point based data with dual phase
        >>> import damask
        >>> import damask.ConfigMaterial as cm
        >>> cm.load_from_Dream3D('20grains16x16x16.dream3D','SyntheticVolumeDataContainer', 'CellData',
        ...                       'EulerAngles','Phases',['Ferrite','Martensite'])

        """
        root_dir = 'DataContainers'
        hdf = h5py.File(fname,'r')

        orientation_path = path.join(root_dir,base_group,data_group,ori_data)
        if hdf[orientation_path].attrs['TupleDimensions'].shape == (3,):
           grain_orientations = np.array(hdf[orientation_path]).reshape(-1,3,order='F')
        else:
           grain_orientations = np.array(hdf[orientation_path])[1:]

        grain_quats = Rotation.from_Euler_angles(grain_orientations).as_quaternion()

        phase_path  = path.join(root_dir,base_group,data_group,phase_id)
        if hdf[phase_path].attrs['TupleDimensions'].shape == (3,):
           grain_phase = np.array(hdf[phase_path]).reshape(-1,order='F')
        else:
           grain_phase = np.array(hdf[phase_path])[1:]

        grain_phase = grain_phase.reshape(len(grain_phase),)
        phase_name_list = [phase_name[i - 1] for i in grain_phase]

        return ConfigMaterial().material_add(phase=phase_name_list, O = grain_quats)                # noqa


    @property
    def is_complete(self):
        """Check for completeness."""
        ok = True
        for top_level in ['homogenization','phase','material']:
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
        """Check for valid content."""
        ok = True

        if 'phase' in self:
            for k,v in self['phase'].items():
                if 'lattice' in v:
                    try:
                        Orientation(lattice=v['lattice'])
                    except KeyError:
                        print(f"Invalid lattice '{v['lattice']}' in phase '{k}'")
                        ok = False

        if 'material' in self:
            for i,m in enumerate(self['material']):
                if 'constituents' in m:
                    v = 0.0
                    for c in m['constituents']:
                        v += float(c['v'])
                        if 'O' in c:
                            try:
                                Rotation.from_quaternion(c['O'])
                            except ValueError:
                                print(f"Invalid orientation '{c['O']}' in material '{i}'")
                                ok = False
                    if not np.isclose(v,1.0):
                        print(f"Total fraction v = {v} ≠ 1 in material '{i}'")
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

        Returns
        -------
        cfg : damask.ConfigMaterial
            Updated material configuration.

        """
        dup = self.copy()
        for i,m in enumerate(dup['material']):
            if ID is not None and i not in ID: continue
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

        Returns
        -------
        cfg : damask.ConfigMaterial
            Updated material configuration.

        """
        dup = self.copy()
        for i,m in enumerate(dup['material']):
            if ID is not None and i not in ID: continue
            try:
                m['homogenization'] = mapping[m['homogenization']]
            except KeyError:
                continue
        return dup


    def material_add(self,**kwargs):
        """
        Add material entries.

        Parameters
        ----------
        **kwargs
            Key-value pairs.

        Returns
        -------
        cfg : damask.ConfigMaterial
            Updated material configuration.

        Examples
        --------
        >>> import numpy as np
        >>> import damask
        >>> m = damask.ConfigMaterial().material_add(phase          = ['Aluminum','Steel'],
        ...                                          O              = damask.Rotation.from_random(2),
        ...                                          homogenization = 'SX')
        >>> m
        material:
          - constituents:
              - O: [0.577764, -0.146299, -0.617669, 0.513010]
                v: 1.0
                phase: Aluminum
            homogenization: SX
          - constituents:
              - O: [0.184176, 0.340305, 0.737247, 0.553840]
                v: 1.0
                phase: Steel
            homogenization: SX
        homogenization: {}
        phase: {}

        >>> m = damask.ConfigMaterial().material_add(phase          = np.array(['Austenite','Martensite']).reshape(1,2),
        ...                                          O              = damask.Rotation.from_random((2,2)),
        ...                                          v              = np.array([0.2,0.8]).reshape(1,2),
        ...                                          homogenization = ['A','B'])
        >>> m
        material:
          - constituents:
              - phase: Austenite
                O: [0.659802978293224, 0.6953785848195171, 0.22426295326327111, -0.17554139512785227]
                v: 0.2
              - phase: Martensite
                O: [0.49356745891301596, 0.2841806579193434, -0.7487679215072818, -0.339085707289975]
                v: 0.8
            homogenization: A
          - constituents:
              - phase: Austenite
                O: [0.26542221365204055, 0.7268854930702071, 0.4474726435701472, -0.44828201137283735]
                v: 0.2
              - phase: Martensite
                O: [0.6545817158479885, -0.08004812803625233, -0.6226561293931374, 0.4212059104577611]
                v: 0.8
            homogenization: B
        homogenization: {}
        phase: {}

        """
        N,n,shaped = 1,1,{}

        for k,v in kwargs.items():
            shaped[k] = np.array(v)
            s = shaped[k].shape[:-1] if k=='O' else shaped[k].shape
            N = max(N,s[0]) if len(s)>0 else N
            n = max(n,s[1]) if len(s)>1 else n

        mat = [{'constituents':[{} for _ in range(n)]} for _ in range(N)]

        if 'v' not in kwargs:
            shaped['v'] = np.broadcast_to(1/n,(N,n))

        for k,v in shaped.items():
            target = (N,n,4) if k=='O' else (N,n)
            obj = np.broadcast_to(v.reshape(util.shapeshifter(v.shape,target,mode='right')),target)
            for i in range(N):
                if k in ['phase','O','v']:
                    for j in range(n):
                        mat[i]['constituents'][j][k] = obj[i,j].item() if isinstance(obj[i,j],np.generic) else obj[i,j]
                else:
                    mat[i][k] = obj[i,0].item() if isinstance(obj[i,0],np.generic) else obj[i,0]

        dup = self.copy()
        dup['material'] = dup['material'] + mat if 'material' in dup else mat

        return dup
