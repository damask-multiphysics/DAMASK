import re
import copy

import pandas as pd
import numpy as np

from . import util

class Table:
    """Manipulate multi-dimensional spreadsheet-like data."""

    def __init__(self,data,shapes,comments=None):
        """
        New spreadsheet.

        Parameters
        ----------
        data : numpy.ndarray or pandas.DataFrame
            Data. Column labels from a pandas.DataFrame will be replaced.
        shapes : dict with str:tuple pairs
            Shapes of the columns. Example 'F':(3,3) for a deformation gradient.
        comments : str or iterable of str, optional
            Additional, human-readable information.

        """
        comments_ = [comments] if isinstance(comments,str) else comments
        self.comments = [] if comments_ is None else [c for c in comments_]
        self.data = pd.DataFrame(data=data)
        self.shapes = { k:(v,) if isinstance(v,(np.int64,np.int32,int)) else v for k,v in shapes.items() }
        self._relabel('uniform')


    def __repr__(self):
        """Brief overview."""
        self._relabel('shapes')
        data_repr = self.data.__repr__()
        self._relabel('uniform')
        return '\n'.join(['# '+c for c in self.comments])+'\n'+data_repr


    def __getitem__(self,item):
        """
        Slice the Table according to item.

        Parameters
        ----------
        item : row and/or column indexer
            Slice to select from Table.

        Returns
        -------
        slice : Table
            Sliced part of the Table.

        Examples
        --------
        >>> import damask
        >>> import numpy as np
        >>> tbl = damask.Table(data=np.arange(12).reshape((4,3)),
        ...                    shapes=dict(colA=(1,),colB=(1,),colC=(1,)))
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
        >>> tbl[1:2,'colB']
           colB
        1     4
        2     7

        """
        item = (item,slice(None,None,None)) if isinstance(item,slice) else \
               item if isinstance(item[0],slice) else \
               (slice(None,None,None),item)
        sliced = self.data.loc[item]
        cols = np.array(sliced.columns if isinstance(sliced,pd.core.frame.DataFrame) else [item[1]])
        _,idx = np.unique(cols,return_index=True)
        return self.__class__(data=sliced,
                              shapes = {k:self.shapes[k] for k in cols[np.sort(idx)]},
                              comments=self.comments)


    def __len__(self):
        """Number of rows."""
        return len(self.data)


    def __copy__(self):
        """Create deep copy."""
        return copy.deepcopy(self)

    copy = __copy__


    def _label(self,what,how):
        """
        Expand labels according to data shape.

        Parameters
        ----------
        what : str or list
            Labels to expand.
        how : str
            Mode of labeling.
            'uniform' ==> v v v
            'shapes'  ==> 3:v v v
            'linear'  ==> 1_v 2_v 3_v

        """
        what = [what] if isinstance(what,str) else what
        labels = []
        for label in what:
            shape = self.shapes[label]
            size = np.prod(shape,dtype=int)
            if   how == 'uniform':
                labels += [label] * size
            elif how == 'shapes':
                labels += [('' if size == 1 or i>0 else f'{util.srepr(shape,"x")}:')+label for i in range(size)]
            elif how == 'linear':
                labels += [('' if size == 1 else f'{i+1}_')+label for i in range(size)]
            else:
                raise KeyError
        return labels


    def _relabel(self,how):
        """
        Modify labeling of data in-place.

        Parameters
        ----------
        how : str
            Mode of labeling.
            'uniform' ==> v v v
            'shapes'  ==> 3:v v v
            'linear'  ==> 1_v 2_v 3_v

        """
        self.data.columns = self._label(self.shapes,how)


    def _add_comment(self,label,shape,info):
        if info is not None:
            specific = f'{label}{" "+str(shape) if np.prod(shape,dtype=int) > 1 else ""}: {info}'
            general  = util.execution_stamp('Table')
            self.comments.append(f'{specific} / {general}')


    def isclose(self,other,rtol=1e-5,atol=1e-8,equal_nan=True):
        """
        Report where values are approximately equal to corresponding ones of other Table.

        Parameters
        ----------
        other : Table
            Table to compare against.
        rtol : float, optional
            Relative tolerance of equality.
        atol : float, optional
            Absolute tolerance of equality.
        equal_nan : bool, optional
            Consider matching NaN values as equal. Defaults to True.

        Returns
        -------
        mask : numpy.ndarray bool
            Mask indicating where corresponding table values are close.

        """
        return np.isclose( self.data.to_numpy(),
                          other.data.to_numpy(),
                          rtol=rtol,
                          atol=atol,
                          equal_nan=equal_nan)


    def allclose(self,other,rtol=1e-5,atol=1e-8,equal_nan=True):
        """
        Test whether all values are approximately equal to corresponding ones of other Table.

        Parameters
        ----------
        other : Table
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
    def load(fname):
        """
        Load from ASCII table file.

        Initial comments are marked by '#', the first non-comment line
        containing the column labels.

        - Vector data column labels are indicated by '1_v, 2_v, ..., n_v'.
        - Tensor data column labels are indicated by '3x3:1_T, 3x3:2_T, ..., 3x3:9_T'.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file for reading.

        """
        try:
            f = open(fname)
        except TypeError:
            f = fname
            f.seek(0)

        comments = []
        line = f.readline().strip()
        while line.startswith('#'):
            comments.append(line.lstrip('#').strip())
            line = f.readline().strip()
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

        return Table(data,shapes,comments)


    @staticmethod
    def load_ang(fname):
        """
        Load from ang file.

        A valid TSL ang file has to have the following columns:

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
            Filename or file for reading.

        """
        try:
            f = open(fname)
        except TypeError:
            f = fname
            f.seek(0)

        content = f.readlines()

        comments = [util.execution_stamp('Table','from_ang')]
        for line in content:
            if line.startswith('#'):
                comments.append(line.split('#',1)[1].strip())
            else:
                break

        data = np.loadtxt(content)

        shapes = {'eu':3, 'pos':2, 'IQ':1, 'CI':1, 'ID':1, 'intensity':1, 'fit':1}
        remainder = data.shape[1]-sum(shapes.values())
        if remainder > 0:                                                       # 3.8 can do: if (remainder := data.shape[1]-sum(shapes.values())) > 0
            shapes['unknown'] = remainder

        return Table(data,shapes,comments)


    @property
    def labels(self):
        return list(self.shapes)


    def get(self,label):
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


    def set(self,label,data,info=None):
        """
        Set column data.

        Parameters
        ----------
        label : str
            Column label.
        data : np.ndarray
            New data.
        info : str, optional
            Human-readable information about the new data.

        Returns
        -------
        table : Table
            Updated table.

        """
        dup = self.copy()
        dup._add_comment(label,data.shape[1:],info)
        m = re.match(r'(.*)\[((\d+,)*(\d+))\]',label)
        if m:
            key = m.group(1)
            idx = np.ravel_multi_index(tuple(map(int,m.group(2).split(","))),
                                       self.shapes[key])
            iloc = dup.data.columns.get_loc(key).tolist().index(True) + idx
            dup.data.iloc[:,iloc] = data
        else:
            dup.data[label]       = data.reshape(dup.data[label].shape)
        return dup


    def add(self,label,data,info=None):
        """
        Add column data.

        Parameters
        ----------
        label : str
            Column label.
        data : np.ndarray
            Modified data.
        info : str, optional
            Human-readable information about the modified data.

        Returns
        -------
        table : Table
            Updated table.

        """
        dup = self.copy()
        dup._add_comment(label,data.shape[1:],info)

        dup.shapes[label] = data.shape[1:] if len(data.shape) > 1 else (1,)
        size = np.prod(data.shape[1:],dtype=int)
        new = pd.DataFrame(data=data.reshape(-1,size),
                           columns=[label]*size,
                          )
        new.index = dup.data.index
        dup.data = pd.concat([dup.data,new],axis=1)
        return dup


    def delete(self,label):
        """
        Delete column data.

        Parameters
        ----------
        label : str
            Column label.

        Returns
        -------
        table : Table
            Updated table.

        """
        dup = self.copy()
        dup.data.drop(columns=label,inplace=True)
        del dup.shapes[label]
        return dup


    def rename(self,old,new,info=None):
        """
        Rename column data.

        Parameters
        ----------
        label_old : str or iterable of str
            Old column label(s).
        label_new : str or iterable of str
            New column label(s).

        Returns
        -------
        table : Table
            Updated table.

        """
        dup = self.copy()
        columns = dict(zip([old] if isinstance(old,str) else old,
                           [new] if isinstance(new,str) else new))
        dup.data.rename(columns=columns,inplace=True)
        dup.comments.append(f'{old} => {new}'+('' if info is None else f': {info}'))
        dup.shapes = {(label if label not in columns else columns[label]):dup.shapes[label] for label in dup.shapes}
        return dup


    def sort_by(self,labels,ascending=True):
        """
        Sort table by values of given labels.

        Parameters
        ----------
        label : str or list
            Column labels for sorting.
        ascending : bool or list, optional
            Set sort order.

        Returns
        -------
        table : Table
            Updated table.

        """
        labels_ = [labels] if isinstance(labels,str) else labels.copy()
        for i,l in enumerate(labels_):
            m = re.match(r'(.*)\[((\d+,)*(\d+))\]',l)
            if m:
                idx = np.ravel_multi_index(tuple(map(int,m.group(2).split(','))),
                                           self.shapes[m.group(1)])
                labels_[i] = f'{1+idx}_{m.group(1)}'

        dup = self.copy()
        dup._relabel('linear')
        dup.data.sort_values(labels_,axis=0,inplace=True,ascending=ascending)
        dup._relabel('uniform')
        dup.comments.append(f'sorted {"ascending" if ascending else "descending"} by {labels}')
        return dup


    def append(self,other):
        """
        Append other table vertically (similar to numpy.vstack).

        Requires matching labels/shapes and order.

        Parameters
        ----------
        other : Table
            Table to append.

        Returns
        -------
        table : Table
            Concatenated table.

        """
        if self.shapes != other.shapes or not self.data.columns.equals(other.data.columns):
            raise KeyError('Labels or shapes or order do not match')
        else:
            dup = self.copy()
            dup.data = dup.data.append(other.data,ignore_index=True)
            return dup


    def join(self,other):
        """
        Append other table horizontally (similar to numpy.hstack).

        Requires matching number of rows and no common labels.

        Parameters
        ----------
        other : Table
            Table to join.

        Returns
        -------
        table : Table
            Joined table.

        """
        if set(self.shapes) & set(other.shapes) or self.data.shape[0] != other.data.shape[0]:
            raise KeyError('Duplicated keys or row count mismatch')
        else:
            dup = self.copy()
            dup.data = dup.data.join(other.data)
            for key in other.shapes:
                dup.shapes[key] = other.shapes[key]
            return dup


    def save(self,fname):
        """
        Save as plain text file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file for writing.

        """
        seen = set()
        labels = []
        for l in [x for x in self.data.columns if not (x in seen or seen.add(x))]:
            if self.shapes[l] == (1,):
                labels.append(f'{l}')
            elif len(self.shapes[l]) == 1:
                labels += [f'{i+1}_{l}' \
                          for i in range(self.shapes[l][0])]
            else:
                labels += [f'{util.srepr(self.shapes[l],"x")}:{i+1}_{l}' \
                          for i in range(np.prod(self.shapes[l]))]

        try:
            fhandle = open(fname,'w',newline='\n')
        except TypeError:
            fhandle = fname

        fhandle.write('\n'.join([f'# {c}' for c in self.comments] + [' '.join(labels)])+'\n')
        self.data.to_csv(fhandle,sep=' ',na_rep='nan',index=False,header=False)
