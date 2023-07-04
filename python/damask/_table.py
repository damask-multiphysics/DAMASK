import re
import copy
from typing import Optional, Union, Tuple, List, Iterable

import pandas as pd
import numpy as np

from ._typehints import FileHandle
from . import util

class Table:
    """Manipulate multi-dimensional spreadsheet-like data."""

    def __init__(self,
                 shapes: dict = {},
                 data: Optional[np.ndarray] = None,
                 comments: Union[None, str, Iterable[str]] = None):
        """
        New spreadsheet.

        Parameters
        ----------
        shapes : dict with str:tuple pairs, optional
            Shapes of the data columns. Mandatory if 'data' is given.
            For instance, 'F':(3,3) for a deformation gradient, or 'r':(1,) for a scalar.
        data : numpy.ndarray or pandas.DataFrame, optional
            Data. Existing column labels of a pandas.DataFrame will be replaced.
        comments : (iterable of) str, optional
            Additional, human-readable information.

        """
        self.comments = [] if comments is None else  \
                        [comments] if isinstance(comments,str) else \
                        [str(c) for c in comments]
        self.shapes = { k:(v,) if isinstance(v,(np.int64,np.int32,int)) else v for k,v in shapes.items() }
        self.data = pd.DataFrame(data=data)
        self._relabel('uniform')


    def __repr__(self) -> str:
        """
        Return repr(self).

        Give short, human-readable summary.

        """
        self._relabel('shapes')
        data_repr = self.data.__repr__()
        self._relabel('uniform')
        return '\n'.join(['# '+c for c in self.comments])+'\n'+data_repr


    def __eq__(self,
               other: object) -> bool:
        """
        Return self==other.

        Test equality of other.

        """
        return NotImplemented if not isinstance(other,Table) else \
               self.shapes == other.shapes and self.data.equals(other.data)


    def __getitem__(self,
                    item: Union[slice, Tuple[slice, ...]]) -> 'Table':
        """
        Return self[item].

        Return slice according to item.

        Parameters
        ----------
        item : row and/or column indexer
            Slice to select from Table.

        Returns
        -------
        slice : damask.Table
            Sliced part of the Table.

        Examples
        --------
        >>> import damask
        >>> import numpy as np
        >>> tbl = damask.Table(shapes=dict(colA=(1,),colB=(1,),colC=(1,)),
        ...                    data=np.arange(12).reshape((4,3)))
        >>> tbl['colA','colB']
           colA  colB
        0     0     1
        1     3     4
        2     6     7
        3     9    10
        >>> tbl[::2,['colB','colA']]
           colB  colA
        0     1     0
        2     7     6
        >>> tbl[[True,False,False,True],'colB']
           colB
        0     1
        3    10

        """
        item_ = (item,slice(None,None,None)) if isinstance(item,(slice,np.ndarray)) else \
                (np.array(item),slice(None,None,None)) if isinstance(item,list) and np.array(item).dtype == np.bool_ else \
                (np.array(item[0]),item[1]) if isinstance(item[0],list) else \
                item if isinstance(item[0],(slice,np.ndarray)) else \
                (slice(None,None,None),item)
        sliced = self.data.loc[item_]
        cols = np.array(sliced.columns if isinstance(sliced,pd.core.frame.DataFrame) else [item_[1]])
        _,idx = np.unique(cols,return_index=True)
        return self.__class__(data=sliced,
                              shapes={k:self.shapes[k] for k in cols[np.sort(idx)]},
                              comments=self.comments)


    def __len__(self) -> int:
        """
        Return len(self).

        Number of rows.

        """
        return len(self.data)


    def __copy__(self) -> 'Table':
        """
        Return deepcopy(self).

        Create deep copy.

        """
        return copy.deepcopy(self)

    copy = __copy__


    def _label(self,
               what: Union[str, List[str]],
               how: str) -> List[str]:
        """
        Expand labels according to data shape.

        Parameters
        ----------
        what : str or list
            Labels to expand.
        how : {'uniform, 'shapes', 'linear'}
            Mode of labeling.
            'uniform' ==> v v v
            'shapes'  ==> 3:v v v
            'linear'  ==> 1_v 2_v 3_v

        """
        what = [what] if isinstance(what,str) else what
        labels = []
        for label in what:
            shape = self.shapes[label]
            size = np.prod(shape,dtype=np.int64)
            if   how == 'uniform':
                labels += [label] * size
            elif how == 'shapes':
                labels += [('' if size == 1 or i>0 else f'{util.srepr(shape,"x")}:')+label for i in range(size)]
            elif how == 'linear':
                labels += [('' if size == 1 else f'{i+1}_')+label for i in range(size)]
            else:
                raise KeyError
        return labels


    def _relabel(self,
                 how: str):
        """
        Modify labeling of data in-place.

        Parameters
        ----------
        how : {'uniform, 'shapes', 'linear'}
            Mode of labeling.
            'uniform' ==> v v v
            'shapes'  ==> 3:v v v
            'linear'  ==> 1_v 2_v 3_v

        """
        self.data.columns = self._label(self.shapes,how)                                            # type: ignore


    def isclose(self,
                other: 'Table',
                rtol: float = 1e-5,
                atol: float = 1e-8,
                equal_nan: bool = True) -> np.ndarray:
        """
        Report where values are approximately equal to corresponding ones of other Table.

        Parameters
        ----------
        other : damask.Table
            Table to compare against.
        rtol : float, optional
            Relative tolerance of equality.
        atol : float, optional
            Absolute tolerance of equality.
        equal_nan : bool, optional
            Consider matching NaN values as equal. Defaults to True.

        Returns
        -------
        mask : numpy.ndarray of bool
            Mask indicating where corresponding table values are close.

        """
        return np.isclose( self.data.to_numpy(),
                          other.data.to_numpy(),
                          rtol=rtol,
                          atol=atol,
                          equal_nan=equal_nan)


    def allclose(self,
                 other: 'Table',
                 rtol: float = 1e-5,
                 atol: float = 1e-8,
                 equal_nan: bool = True) -> bool:
        """
        Test whether all values are approximately equal to corresponding ones of other Table.

        Parameters
        ----------
        other : damask.Table
            Table to compare against.
        rtol : float, optional
            Relative tolerance of equality.
        atol : float, optional
            Absolute tolerance of equality.
        equal_nan : bool, optional
            Consider matching NaN values as equal. Defaults to True.

        Returns
        -------
        answer : bool
            Whether corresponding values are close between both tables.

        """
        return np.allclose( self.data.to_numpy(),
                           other.data.to_numpy(),
                           rtol=rtol,
                           atol=atol,
                           equal_nan=equal_nan)


    @staticmethod
    def load(fname: FileHandle) -> 'Table':
        """
        Load from ASCII table file.

        Initial comments are marked by '#'.
        The first non-comment line contains the column labels.

        - Vector data column labels are indicated by '1_v, 2_v, ..., n_v'.
        - Tensor data column labels are indicated by '3x3:1_T, 3x3:2_T, ..., 3x3:9_T'.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file to read.

        Returns
        -------
        loaded : damask.Table
            Table data from file.

        """
        f = util.open_text(fname)
        f.seek(0)

        comments = []
        while (line := f.readline().strip()).startswith('#'):
            comments.append(line.lstrip('#').strip())
        labels = line.split()

        shapes = {}
        for label in labels:
            tensor_column = re.search(r'[0-9,x]*?:[0-9]*?_',label)
            if tensor_column:
                my_shape = tensor_column.group().split(':',1)[0].split('x')
                shapes[label.split('_',1)[1]] = tuple([int(d) for d in my_shape])
            else:
                vector_column = re.match(r'[0-9]*?_',label)
                if vector_column:
                    shapes[label.split('_',1)[1]] = (int(label.split('_',1)[0]),)
                else:
                    shapes[label] = (1,)

        data = pd.read_csv(f,names=list(range(len(labels))),sep=r'\s+')

        return Table(shapes,data,comments)


    @staticmethod
    def load_ang(fname: FileHandle,
                 shapes = {'eu':3,
                           'pos':2,
                           'IQ':1,
                           'CI':1,
                           'ID':1,
                           'intensity':1,
                           'fit':1}) -> 'Table':
        """
        Load from ANG file.

        Regular ANG files feature the following columns:

        - Euler angles (Bunge notation) in radians, 3 floats, label 'eu'.
        - Spatial position in meters, 2 floats, label 'pos'.
        - Image quality, 1 float, label 'IQ'.
        - Confidence index, 1 float, label 'CI'.
        - Phase ID, 1 int, label 'ID'.
        - SEM signal, 1 float, label 'intensity'.
        - Fit, 1 float, label 'fit'.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file to read.
        shapes : dict with str:int pairs, optional
            Column labels and their width.
            Defaults to standard TSL ANG format.

        Returns
        -------
        loaded : damask.Table
            Table data from file.

        """
        f = util.open_text(fname)
        f.seek(0)

        content = f.readlines()

        comments = [util.execution_stamp('Table','from_ang')]
        for line in content:
            if line.startswith('#'):
                comments.append(line.split('#',1)[1].strip())
            else:
                break

        data = np.loadtxt(content)

        if (remainder := data.shape[1]-sum(shapes.values())) > 0:
            shapes['unknown'] = remainder

        return Table(shapes,data,comments)


    @property
    def labels(self) -> List[str]:
        return list(self.shapes)


    def get(self,
            label: str) -> np.ndarray:
        """
        Get column data.

        Parameters
        ----------
        label : str
            Column label.

        Returns
        -------
        data : numpy.ndarray
            Array of column data.

        """
        data = self.data[label].to_numpy().reshape((-1,)+self.shapes[label])

        return data.astype(type(data.flatten()[0]))


    def set(self,
            label: str,
            data: np.ndarray,
            info: Optional[str] = None) -> 'Table':
        """
        Add new or replace existing column data.

        Parameters
        ----------
        label : str
            Column label.
        data : numpy.ndarray
            Column data. First dimension needs to match number of rows.
        info : str, optional
            Human-readable information about the data.

        Returns
        -------
        updated : damask.Table
            Updated table.

        """
        def add_comment(label: str, shape: Tuple[int, ...],info: str) -> List[str]:
            specific = f'{label}{" "+str(shape) if np.prod(shape,dtype=np.int64) > 1 else ""}: {info}'
            general  = util.execution_stamp('Table')
            return [f'{specific} / {general}']

        dup = self.copy()
        if info is not None: self.comments += add_comment(label,data.shape[1:],info)

        key = m.group(1) if (m := re.match(r'(.*)\[((\d+,)*(\d+))\]',label)) else label

        if key in dup.shapes:

            if m:
                idx = np.ravel_multi_index(tuple(map(int,m.group(2).split(","))),
                                           self.shapes[key])
                iloc = dup.data.columns.get_loc(key).tolist().index(True) + idx
                dup.data.iloc[:,iloc] = data
            else:
                dup.data[label]       = data.reshape(dup.data[label].shape)

        else:

            dup.shapes[label] = data.shape[1:] if len(data.shape) > 1 else (1,)
            size = np.prod(data.shape[1:],dtype=np.int64)
            new = pd.DataFrame(data=data.reshape(-1,size),
                               columns=[label]*size,
                              )
            new.index = new.index if dup.data.index.empty else dup.data.index
            dup.data = pd.concat([dup.data,new],axis=1)

        return dup


    def delete(self,
               label: str) -> 'Table':
        """
        Delete column data.

        Parameters
        ----------
        label : str
            Column label.

        Returns
        -------
        updated : damask.Table
            Updated table.

        """
        dup = self.copy()
        dup.data.drop(columns=label,inplace=True)
        del dup.shapes[label]
        return dup


    def rename(self,
               old: Union[str, Iterable[str]],
               new: Union[str, Iterable[str]],
               info: Optional[str] = None) -> 'Table':
        """
        Rename column data.

        Parameters
        ----------
        label_old : (iterable of) str
            Old column labels.
        label_new : (iterable of) str
            New column labels.

        Returns
        -------
        updated : damask.Table
            Updated table.

        """
        dup = self.copy()
        columns = dict(zip([old] if isinstance(old,str) else old,
                           [new] if isinstance(new,str) else new))
        dup.data.rename(columns=columns,inplace=True)
        dup.comments.append(f'{old} => {new}'+('' if info is None else f': {info}'))
        dup.shapes = {(label if label not in columns else columns[label]):dup.shapes[label] for label in dup.shapes}
        return dup


    def sort_by(self,
                labels: Union[str, List[str]],
                ascending: Union[bool, List[bool]] = True) -> 'Table':
        """
        Sort table by data of given columns.

        Parameters
        ----------
        label : str or list
            Column labels for sorting.
        ascending : bool or list, optional
            Set sort order. Defaults to True.

        Returns
        -------
        updated : damask.Table
            Updated table.

        """
        labels_ = [labels] if isinstance(labels,str) else labels.copy()
        for i,l in enumerate(labels_):
            if m := re.match(r'(.*)\[((\d+,)*(\d+))\]',l):
                idx = np.ravel_multi_index(tuple(map(int,m.group(2).split(','))),
                                           self.shapes[m.group(1)])
                labels_[i] = f'{1+idx}_{m.group(1)}'

        dup = self.copy()
        dup._relabel('linear')
        dup.data.sort_values(labels_,axis=0,inplace=True,ascending=ascending)
        dup._relabel('uniform')
        dup.comments.append(f'sorted {"ascending" if ascending else "descending"} by {labels}')
        return dup


    def append(self,
               other: 'Table') -> 'Table':
        """
        Append other table vertically (similar to numpy.vstack).

        Requires matching labels/shapes and order.

        Parameters
        ----------
        other : damask.Table
            Table to append.

        Returns
        -------
        updated : damask.Table
            Updated table.

        """
        if self.shapes != other.shapes or not self.data.columns.equals(other.data.columns):
            raise KeyError('mismatch of shapes or labels or their order')

        dup = self.copy()
        dup.data = pd.concat([dup.data,other.data],ignore_index=True)
        return dup


    def join(self,
             other: 'Table') -> 'Table':
        """
        Append other table horizontally (similar to numpy.hstack).

        Requires matching number of rows and no common labels.

        Parameters
        ----------
        other : damask.Table
            Table to join.

        Returns
        -------
        updated : damask.Table
            Updated table.

        """
        if set(self.shapes) & set(other.shapes) or self.data.shape[0] != other.data.shape[0]:
            raise KeyError('duplicated keys or row count mismatch')

        dup = self.copy()
        dup.data = dup.data.join(other.data)
        for key in other.shapes:
            dup.shapes[key] = other.shapes[key]
        return dup


    def save(self,
             fname: FileHandle,
             with_labels: bool = True):
        """
        Save as plain text file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file to write.
        with_labels : bool, optional
            Write column labels. Defaults to True.

        """
        labels = []
        if with_labels:
            for l in list(dict.fromkeys(self.data.columns)):
                if self.shapes[l] == (1,):
                    labels.append(f'{l}')
                elif len(self.shapes[l]) == 1:
                    labels += [f'{i+1}_{l}' \
                              for i in range(self.shapes[l][0])]
                else:
                    labels += [f'{util.srepr(self.shapes[l],"x")}:{i+1}_{l}' \
                              for i in range(np.prod(self.shapes[l]))]

        f = util.open_text(fname,'w')

        f.write('\n'.join([f'# {c}' for c in self.comments] + [' '.join(labels)])+('\n' if labels else ''))
        try:                                                                                        # backward compatibility
            self.data.to_csv(f,sep=' ',na_rep='nan',index=False,header=False,lineterminator='\n')
        except TypeError:
            self.data.to_csv(f,sep=' ',na_rep='nan',index=False,header=False,line_terminator='\n')
