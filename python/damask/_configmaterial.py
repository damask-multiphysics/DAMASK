from typing import Optional, Union, Sequence, Dict, Any, List

import numpy as np
import h5py
import logging

from . import YAML
from . import Rotation
from . import Orientation
from . import util
from . import tensor
from . import Table
from ._typehints import FloatSequence, StrSequence


logger = logging.getLogger(__name__)

class ConfigMaterial(YAML):
    """
    Material configuration.

    Manipulate material configurations for storage in YAML format.
    A complete material configuration file has the entries 'material',
    'phase', and 'homogenization'.
    """

    def __init__(self,
                 config: Optional[Union[str,Dict[str,Any]]] = None,*,
                 homogenization: Optional[Dict[str,Dict]] = None,
                 phase: Optional[Dict[str,Dict]] = None,
                 material: Optional[List[Dict[str,Any]]] = None):
        """
        New material configuration.

        Parameters
        ----------
        config : dict or str, optional
            Material configuration. String needs to be valid YAML.
        homogenization : dict, optional
            Homogenization configuration.
            Defaults to an empty dict if 'config' is not given.
        phase : dict, optional
            Phase configuration.
            Defaults to an empty dict if 'config' is not given.
        material : dict, optional
            Materialpoint configuration.
            Defaults to an empty list if 'config' is not given.
        """
        kwargs: Dict[str,Union[Dict[str,Dict],List[Dict[str,Any]]]] = {}
        for arg,value in zip(['homogenization','phase','material'],[homogenization,phase,material]):
            if value is None and config is None:
                kwargs[arg] = [] if arg == 'material' else {}
            elif value is not None:
                kwargs[arg] = value

        super().__init__(config,**kwargs)


    @staticmethod
    def load_DREAM3D(fname: str,
                     grain_data: Optional[str] = None,
                     cell_data: Optional[str] = None,
                     cell_ensemble_data: str = 'CellEnsembleData',
                     phases: str = 'Phases',
                     Euler_angles: str = 'EulerAngles',
                     phase_names: str = 'PhaseName',
                     base_group: Optional[str] = None) -> 'ConfigMaterial':
        """
        Load DREAM.3D (HDF5) file.

        Data in DREAM.3D files can be stored per cell ('CellData')
        and/or per grain ('Grain Data'). Per default, i.e. if
        'grain_data' is None, cell-wise data is assumed.

        Parameters
        ----------
        fname : str or pathlib.Path
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

        Returns
        -------
        loaded : damask.ConfigMaterial
            Material configuration from file.

        Notes
        -----
        A grain-wise material configuration is based on segmented data from
        the DREAM.3D file. This data is typically available when the microstructure
        was synthetically created. In cell-wise representations, cells having the
        same orientation and phase are grouped. Since synthetically created
        microstructures have typically no in-grain scatter, cell-wise grids
        can appear to be segmented.

        damask.GeomGrid.load_DREAM3D creates the corresponding grid-based
        geometry definition. Since the numbering of materials in cell-wise
        and grain-wise grids is different, it is imperative to use the same
        mode for both load_DREAM3D functions. That means, if the "grain_data"
        argument is used for this function, the correct grid configuration
        is only obtained if the "feature_IDs" argument is used when calling
        damask.GeomGrid.load_DREAM3D.

        Homogenization and phase entries are emtpy and need to be
        defined separately.
        """
        with h5py.File(fname, 'r') as f:
            b = util.DREAM3D_base_group(f) if base_group is None else base_group
            c = util.DREAM3D_cell_data_group(f) if cell_data is None else cell_data

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

        base_config = ConfigMaterial({'phase':{k if isinstance(k,int) else str(k): None for k in np.unique(phase)},
                                      'homogenization':{'direct':{'N_constituents':1}}})
        constituent = {k:np.atleast_1d(v[idx].squeeze()) for k,v in zip(['O','phase'],[O,phase])}

        return base_config.material_add(**constituent,homogenization='direct')


    @staticmethod
    def from_table(table: Table,*,
                   homogenization: Optional[Union[str,StrSequence]] = None,
                   phase: Optional[Union[str,StrSequence]] = None,
                   v: Optional[Union[str,FloatSequence]] = None,
                   O: Optional[Union[str,FloatSequence]] = None,
                   V_e: Optional[Union[str,FloatSequence]] = None) -> 'ConfigMaterial':
        """
        Generate from an ASCII table.

        Parameters
        ----------
        table : damask.Table
            Table that contains material information.
        homogenization : (array-like) of str, optional
            Homogenization label.
        phase : (array-like) of str, optional
            Phase label (per constituent).
        v : (array-like) of float or str, optional
            Constituent volume fraction (per constituent).
            Defaults to 1/N_constituent.
        O : (array-like) of damask.Rotation or np.array/list of shape (4) or str, optional
            Orientation as unit quaternion (per constituent).
        V_e : (array-like) of np.array/list of shape (3,3) or str, optional
            Left elastic stretch (per constituent).

        Returns
        -------
        new : damask.ConfigMaterial
            Material configuration from values in table.

        Notes
        -----
        If the value of an argument is a string that is a column label,
        data from the table is used to fill the corresponding entry in
        the material configuration. Otherwise, the value is used directly.

        First index of array-like values that are defined per constituent
        runs over materials, whereas second index runs over constituents.

        Examples
        --------
        >>> import damask
        >>> from damask import ConfigMaterial as cm
        >>> t = damask.Table.load('small.txt')
        >>> t
            3:pos  pos  pos  4:qu   qu   qu   qu     phase homog
         0      0    0    0   1.0  0.0  0.0  0.0  Aluminum    SX
         1      1    0    0   0.0  1.0  0.0  0.0     Steel    SX
         2      1    1    0   0.0  1.0  0.0  0.0     Steel    SX

        >>> cm.from_table(t,O='qu',phase='phase',homogenization='homog')
        homogenization: {SX: null}
        phase: {Aluminum: null, Steel: null}
        material:
        - constituents:
          - phase: Aluminum
            O: [1.0, 0.0, 0.0, 0.0]
            v: 1.0
          homogenization: SX
        - constituents:
          - phase: Steel
            O: [0.0, 1.0, 0.0, 0.0]
            v: 1.0
          homogenization: SX

        >>> cm.from_table(t,O='qu',phase='phase',homogenization='single_crystal')
        homogenization: {single_crystal: null}
        phase: {Aluminum: null, Steel: null}
        material:
        - constituents:
          - phase: Aluminum
            O: [1.0, 0.0, 0.0, 0.0]
            v: 1.0
          homogenization: single_crystal
        - constituents:
          - phase: Steel
            O: [0.0, 1.0, 0.0, 0.0]
            v: 1.0
          homogenization: single_crystal
        """
        kwargs = {}
        for arg,val in zip(['homogenization','phase','v','O','V_e'],[homogenization,phase,v,O,V_e]):
            if val is not None:
                kwargs[arg] = table.get(val) if val in table.labels else np.atleast_2d([val]*len(table)).T # type: ignore[arg-type]

        _,idx = np.unique(np.hstack(list(kwargs.values())),return_index=True,axis=0)
        idx = np.sort(idx)
        kwargs = {k:np.atleast_1d(v[idx].squeeze()) for k,v in kwargs.items()}

        return ConfigMaterial().material_add(**kwargs)


    @property
    def is_complete(self) -> bool:
        """
        Check for completeness.

        Only the general file layout is considered.
        This check does not consider whether specific parameters for
        a particular phase/homogenization model are missing.

        Returns
        -------
        complete : bool
            Whether the material.yaml definition is complete.
        """
        def LabeledList(label,items):
            return f'{label.capitalize()}{"s" if len(items)>1 else ""} {util.srepr(items,",",quote=True)}'

        ok = True
        msg = []
        all = set(['homogenization','phase','material'])
        miss = set([item for item in all if item not in self])
        empty = set([item for item in all-miss if self[item] is None])

        if miss:
            msg.append(f'{LabeledList("top-level",miss)} missing')
            ok = False
        if empty:
            msg.append(f'{LabeledList("top-level",empty)} empty')

        if ok:
            ok &= len(self['material']) > 0
            if len(self['material']) < 1: msg.append('No materials defined')

            homogenization = set()
            phase          = set()
            for i,v in enumerate(self['material']):
                if 'homogenization' in v:
                    homogenization.add(v['homogenization'])
                else:
                    msg.append(f'No homogenization specified for material {i}')
                    ok = False

                if 'constituents' in v:
                    for ii,vv in enumerate(v['constituents']):
                        if 'O' not in vv:
                            msg.append(f'No orientation specified for constituent {ii} of material {i}')
                            ok = False
                        if 'phase' in vv:
                            phase.add(vv['phase'])
                        else:
                            msg.append(f'No phase specified for constituent {ii} of material {i}')
                            ok = False

            for v,other in {'phase':phase,
                            'homogenization':homogenization}.items():
                me = set([] if v in empty else self[v])
                if _miss := other - me:
                    msg.append(f'{LabeledList(v,_miss)} missing')
                    ok = False
                if len(_empty := [item for item in me if self[v][item] is None]) > 0:
                    msg.append(f'{LabeledList(v,_empty)} undefined')
                    ok = False

        logger.info(util.srepr(msg))
        return ok


    @property
    def is_valid(self) -> bool:
        """
        Check for valid content.

        Only the generic file content is considered.
        This check does not consider whether parameters for a
        particular phase/homogenization mode are out of bounds.

        Returns
        -------
        valid : bool
            Whether the material.yaml definition is valid.
        """
        ok = True

        if 'phase' in self:
            for k,v in self['phase'].items():
                if v is not None and 'lattice' in v:
                    try:
                        Orientation(lattice=v['lattice'])
                    except KeyError:
                        logger.warning(f"Invalid lattice '{v['lattice']}' in phase '{k}'")
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
                                logger.warning(f"Invalid orientation '{c['O']}' in material '{i}'")
                                ok = False
                    if not np.isclose(v,1.0):
                        logger.warning(f"Total fraction v = {v} ≠ 1 in material '{i}'")
                        ok = False

        return ok


    def material_rename_phase(self,
                              mapping: Dict[str, str],
                              ID: Optional[Sequence[int]] = None,
                              constituent: Optional[Sequence[int]] = None) -> 'ConfigMaterial':
        """
        Change phase name in material.

        Parameters
        ----------
        mapping : dict
            Mapping from old name to new name.
        ID : list of ints, optional
            Limit renaming to selected material IDs.
        constituent : list of ints, optional
            Limit renaming to selected constituents.

        Returns
        -------
        updated : damask.ConfigMaterial
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


    def material_rename_homogenization(self,
                                       mapping: Dict[str, str],
                                       ID: Optional[Sequence[int]] = None) -> 'ConfigMaterial':
        """
        Change homogenization name in material.

        Parameters
        ----------
        mapping : dict
            Mapping from old name to new name.
        ID : list of ints, optional
            Limit renaming to selected homogenization IDs.

        Returns
        -------
        updated : damask.ConfigMaterial
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


    def material_add(self,*,
                     homogenization: Optional[Union[str,StrSequence]] = None,
                     phase: Optional[Union[str,StrSequence]] = None,
                     v: Optional[Union[float,FloatSequence]] = None,
                     O: Optional[Union[float,FloatSequence]] = None,
                     V_e: Optional[Union[float,FloatSequence]] = None) -> 'ConfigMaterial':
        """
        Add material entries.

        Parameters
        ----------
        homogenization : (array-like) of str, optional
            Homogenization label.
        phase : (array-like) of str, optional
            Phase label (per constituent).
        v : (array-like) of float, optional
            Constituent volume fraction (per constituent).
            Defaults to 1/N_constituent.
        O : (array-like) of damask.Rotation or np.array/list of shape (4), optional
            Orientation as unit quaternion (per constituent).
        V_e : (array-like) of np.array/list of shape (3,3), optional
            Left elastic stretch (per constituent).

        Returns
        -------
        updated : damask.ConfigMaterial
            Updated material configuration.

        Notes
        -----
        First index of array-like values that are defined per constituent
        runs over materials, whereas second index runs over constituents.

        Examples
        --------
        Create two grains of ferrite and one grain of martensite, each with random orientation:

        >>> import damask
        >>> m = damask.ConfigMaterial()
        >>> m = m.material_add(phase = ['Ferrite','Martensite','Ferrite'],
        ...                    O = damask.Rotation.from_random(3,rng_seed=20191102),
        ...                    homogenization = 'SX')
        >>> m
        homogenization: {SX: null}
        phase: {Ferrite: null, Martensite: null}
        material:
          - constituents:
              - phase: Ferrite
                O: [0.0047, -0.9582, 0.1084, 0.2645]
                v: 1.0
            homogenization: SX
          - constituents:
              - phase: Martensite
                O: [0.9147, -0.1907, 0.2901, -0.2068]
                v: 1.0
            homogenization: SX
          - constituents:
              - phase: Ferrite
                O: [0.1068, -0.4427, 0.1369, 0.8797]
                v: 1.0
            homogenization: SX

        Create five materials that each approximate a duplex stainless steel microstructure
        with three austenite and one relatively bigger ferrite grain of random orientation each:

        >>> import numpy as np
        >>> import damask
        >>> m = damask.ConfigMaterial()
        >>> N_materials = 5
        >>> m = m.material_add(phase = np.array([['Austenite']*3+['Ferrite']]),
        ...                    O = damask.Rotation.from_random((N_materials,4),rng_seed=20191102),
        ...                    v = np.array([[0.2]*3+[0.4]]),
        ...                    homogenization = 'Taylor')
        >>> m
        homogenization: {Taylor: null}
        phase: {Austenite: null, Ferrite: null}
        material:
        - constituents:
          - phase: Austenite
            v: 0.2
            O: [0.004702411137213036, -0.9582446864633862, 0.1084379916089085, 0.2645490694937509]
          - phase: Austenite
            v: 0.2
            O: [0.9147097460704486, -0.19068436891182194, 0.29014401444532145, -0.20678975501215882]
          - phase: Austenite
            v: 0.2
            O: [0.10677819003833185, -0.4427133706883004, 0.13690394495734726, 0.879693468999888]
          - phase: Ferrite
            v: 0.4
            O: [0.8664338002923555, 0.04448357787828491, -0.4945927532088464, 0.05188149461403649]
          homogenization: Taylor
        - constituents:
          - phase: Austenite
            v: 0.2
            O: [0.5621873738314133, 0.0028841916095125584, -0.817023371343172, -0.1281009321680984]
          - phase: Austenite
            v: 0.2
            O: [0.1566777437467901, -0.8117282158019414, 0.5096142534839398, 0.23841707348975383]
          - phase: Austenite
            v: 0.2
            O: [0.3559036203819333, 0.1946923701552408, 0.058744995087853975, -0.9121274689178566]
          - phase: Ferrite
            v: 0.4
            O: [0.467387781713959, -0.35644325887489176, 0.8031986430613528, 0.09679258489963502]
          homogenization: Taylor
        - constituents:
          - phase: Austenite
            v: 0.2
            O: [0.4399087544327661, 0.12802483830067418, -0.8257167208737983, 0.32906203886337354]
          - phase: Austenite
            v: 0.2
            O: [0.12410381094181624, -0.5125024631828828, -0.8493860709598213, 0.021972068647108236]
          - phase: Austenite
            v: 0.2
            O: [0.03909373022192218, 0.4596226773046959, 0.42809626138739537, 0.7771436583738773]
          - phase: Ferrite
            v: 0.4
            O: [0.737821660605232, 0.38809925187040367, -0.012129167758963711, 0.5521331824196455]
          homogenization: Taylor
        - constituents:
          - phase: Austenite
            v: 0.2
            O: [0.4924738838478857, -0.0534798919571679, -0.6570981342247908, 0.5681825559468784]
          - phase: Austenite
            v: 0.2
            O: [0.13073521303792138, 0.2534173177988532, -0.9582490914178947, -0.021133998872519554]
          - phase: Austenite
            v: 0.2
            O: [0.1633346595899539, 0.6775968809652247, -0.07127256805012916, -0.71351557581203]
          - phase: Ferrite
            v: 0.4
            O: [0.7658044627436773, -0.5327872540278646, 0.1102330397070761, 0.34282640467772235]
          homogenization: Taylor
        - constituents:
          - phase: Austenite
            v: 0.2
            O: [0.25814496892598815, -0.6159961898524933, -0.5080223627084379, 0.543896265930874]
          - phase: Austenite
            v: 0.2
            O: [0.8497433829153472, 0.4264182767672584, 0.05570674517418605, -0.3049596612218108]
          - phase: Austenite
            v: 0.2
            O: [0.5146112784760113, 0.529467219604771, 0.661078636611197, 0.13347183839881469]
          - phase: Ferrite
            v: 0.4
            O: [0.18430893147208752, 0.012407731059331692, -0.5551804816056372, -0.8109567798802285]
          homogenization: Taylor
        """
        dim = {'O':(4,),'V_e':(3,3,)}
        ex = dict((keyword, -len(val)) for keyword,val in dim.items())

        N_materials,N_constituents = 1,1
        shape = {}
        for arg,val in zip(['homogenization','phase','v','O','V_e'],[homogenization,phase,v,O,V_e]):
            if val is None: continue
            shape[arg] = np.array(val)
            s = shape[arg].shape[:ex.get(arg,None)]
            N_materials = max(N_materials,s[0]) if len(s)>0 else N_materials
            N_constituents = max(N_constituents,s[1]) if len(s)>1 else N_constituents

        shape['v'] = np.array(shape.get('v',1./N_constituents),float)

        mat: Sequence[dict] = [{'constituents':[{} for _ in range(N_constituents)]} for _ in range(N_materials)]

        for k,v in shape.items():
            target = (N_materials,N_constituents) + dim.get(k,())
            broadcasted = np.broadcast_to(np.array(v).reshape(util.shapeshifter(np.array(v).shape,target,'right')),target)
            if k == 'v':
                if np.min(broadcasted) < 0 or np.max(broadcasted) > 1:
                    raise ValueError('volume fraction "v" out of range')
                if len(np.atleast_1d(broadcasted)) > 1:
                    total = np.sum(broadcasted,axis=-1)
                    if np.min(total) < 0 or np.max(total) > 1:
                        raise ValueError('volume fraction "v" out of range')
            if k == 'O' and not np.allclose(1.0,np.linalg.norm(broadcasted,axis=-1)):
                raise ValueError('orientation "O" is not a unit quaterion')
            elif k == 'V_e' and not np.allclose(broadcasted,tensor.symmetric(broadcasted)):
                raise ValueError('elastic stretch "V_e" is not symmetric')
            for i in range(N_materials):
                if k == 'homogenization':
                    mat[i][k] = broadcasted[i,0]
                else:
                    for j in range(N_constituents):
                        mat[i]['constituents'][j][k] = broadcasted[i,j]

        dup = self.copy()
        dup['material'] = dup['material'] + mat if 'material' in dup else mat

        for what in [item for item in ['phase','homogenization'] if item in shape]:
            for k in np.unique(shape[what]):
                if k not in dup[what]: dup[what][str(k)] = None

        return dup
