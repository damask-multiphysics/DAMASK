import numpy as np
import h5py

from . import Config
from . import Rotation
from . import Orientation
from . import util

class ConfigMaterial(Config):
    """
    Material configuration.

    Manipulate material configurations for storage in YAML format.
    A complete material configuration file has the entries 'material',
    'phase', and 'homogenization'. For use in DAMASK, it needs to be
    stored as 'material.yaml'.
    
    """

    def __init__(self,d=None):
        """
        New material configuration.

        Parameters
        ----------
        d : dictionary, optional
            Initial content. Defaults to None, in which case empty entries for
            material, homogenization, and phase are created.

        """
        super().__init__({'material': [], 'homogenization': {}, 'phase': {}} if d is None else d)


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
        Generate from an ASCII table.

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
        kwargs_ = {k:table.get(v) for k,v in kwargs.items()}

        _,idx = np.unique(np.hstack(list(kwargs_.values())),return_index=True,axis=0)
        idx = np.sort(idx)
        kwargs_ = {k:np.atleast_1d(v[idx].squeeze()) for k,v in kwargs_.items()}

        return ConfigMaterial().material_add(**kwargs_)


    @staticmethod
    def load_DREAM3D(fname,
                     grain_data=None,cell_data=None,cell_ensemble_data='CellEnsembleData',
                     phases='Phases',Euler_angles='EulerAngles',phase_names='PhaseName',
                     base_group=None):
        """
        Load DREAM.3D (HDF5) file.

        Data in DREAM.3D files can be stored per cell ('CellData')
        and/or per grain ('Grain Data'). Per default, cell-wise data
        is assumed.

        damask.Grid.load_DREAM3D allows to get the corresponding geometry
        for the grid solver.

        Parameters
        ----------
        fname : str
            Filename of the DREAM.3D (HDF5) file.
        grain_data : str
            Name of the group (folder) containing grain-wise data. Defaults
            to None, in which case cell-wise data is used.
        cell_data : str
            Name of the group (folder) containing cell-wise data. Defaults to
            None in wich case it is automatically detected.
        cell_ensemble_data : str
            Name of the group (folder) containing data of cell ensembles. This
            group is used to inquire the name of the phases. Phases will get
            numeric IDs if this group is not found. Defaults to 'CellEnsembleData'.
        phases : str
            Name of the dataset containing the phase ID (cell-wise or grain-wise).
            Defaults to 'Phases'.
        Euler_angles : str
            Name of the dataset containing the crystallographic orientation as
            Euler angles in radians (cell-wise or grain-wise). Defaults to 'EulerAngles'.
        phase_names : str
            Name of the dataset containing the phase names. Phases will get
            numeric IDs if this dataset is not found. Defaults to 'PhaseName'.
        base_group : str
            Path to the group (folder) that contains geometry (_SIMPL_GEOMETRY),
            and grain- or cell-wise data. Defaults to None, in which case
            it is set as the path that contains _SIMPL_GEOMETRY/SPACING.

        """
        b = util.DREAM3D_base_group(fname) if base_group is None else base_group
        c = util.DREAM3D_cell_data_group(fname) if cell_data is None else cell_data
        f = h5py.File(fname,'r')

        if grain_data is None:
            phase = f['/'.join([b,c,phases])][()].flatten()
            O = Rotation.from_Euler_angles(f['/'.join([b,c,Euler_angles])]).as_quaternion().reshape(-1,4) # noqa
            _,idx = np.unique(np.hstack([O,phase.reshape(-1,1)]),return_index=True,axis=0)
            idx = np.sort(idx)
        else:
            phase = f['/'.join([b,grain_data,phases])][()]
            O = Rotation.from_Euler_angles(f['/'.join([b,grain_data,Euler_angles])]).as_quaternion() # noqa
            idx = np.arange(phase.size)

        if cell_ensemble_data is not None and phase_names is not None:
            try:
                names = np.array([s.decode() for s in f['/'.join([b,cell_ensemble_data,phase_names])]])
                phase = names[phase]
            except KeyError:
                pass


        base_config = ConfigMaterial({'phase':{k if isinstance(k,int) else str(k):'t.b.d.' for k in np.unique(phase)},
                                      'homogenization':{'direct':{'N_constituents':1}}})
        constituent = {k:np.atleast_1d(v[idx].squeeze()) for k,v in zip(['O','phase'],[O,phase])}

        return base_config.material_add(**constituent,homogenization='direct')


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
