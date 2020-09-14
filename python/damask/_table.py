import re
import copy

import pandas as pd
import numpy as np

from . import util

class Table:
    """Store spreadsheet-like data."""

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
        self.shapes = { k:(v,) if isinstance(v,(np.int,int)) else v for k,v in shapes.items() }
        self._label_condensed()


    def __copy__(self):
        """Copy Table."""
        return copy.deepcopy(self)

    def copy(self):
        """Copy Table."""
        return self.__copy__()


    def _label_flat(self):
        """Label data individually, e.g. v v v ==> 1_v 2_v 3_v."""
        labels = []
        for label,shape in self.shapes.items():
            size = int(np.prod(shape))
            labels += [('' if size == 1 else f'{i+1}_')+label for i in range(size)]
        self.data.columns = labels


    def _label_condensed(self):
        """Label data condensed, e.g. 1_v 2_v 3_v ==> v v v."""
        labels = []
        for label,shape in self.shapes.items():
            labels += [label] * int(np.prod(shape))
        self.data.columns = labels


    def _add_comment(self,label,shape,info):
        if info is not None:
            specific = f'{label}{" "+str(shape) if np.prod(shape,dtype=int) > 1 else ""}: {info}'
            general  = util.execution_stamp('Table')
            self.comments.append(f'{specific} / {general}')


    @staticmethod
    def from_ASCII(fname):
        """
        Create table from ASCII file.

        The first line can indicate the number of subsequent header lines as 'n header',
        alternatively first line is the header and comments are marked by '#' ('new style').
        Vector data column labels are indicated by '1_v, 2_v, ..., n_v'.
        Tensor data column labels are indicated by '3x3:1_T, 3x3:2_T, ..., 3x3:9_T'.

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

        try:
            N_comment_lines,keyword = f.readline().strip().split(maxsplit=1)
            if keyword != 'header':
                raise ValueError
            else:
                comments = [f.readline().strip() for i in range(1,int(N_comment_lines))]
                labels   = f.readline().split()
        except ValueError:
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
    def from_ang(fname):
        """
        Create table from TSL ang file.

        A valid TSL ang file needs to contains the following columns:
        * Euler angles (Bunge notation) in radians, 3 floats, label 'eu'.
        * Spatial position in meters, 2 floats, label 'pos'.
        * Image quality, 1 float, label 'IQ'.
        * Confidence index, 1 float, label 'CI'.
        * Phase ID, 1 int, label 'ID'.
        * SEM signal, 1 float, label 'intensity'.
        * Fit, 1 float, label 'fit'.

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
                comments.append(line.strip())
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
        return list(self.shapes.keys())


    def get(self,label):
        """
        Get column data.

        Parameters
        ----------
        label : str
            Column label.

        """
        if re.match(r'[0-9]*?_',label):
            idx,key = label.split('_',1)
            data = self.data[key].to_numpy()[:,int(idx)-1].reshape(-1,1)
        else:
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

        """
        dup = self.copy()
        dup._add_comment(label,data.shape[1:],info)

        if re.match(r'[0-9]*?_',label):
            idx,key = label.split('_',1)
            iloc = dup.data.columns.get_loc(key).tolist().index(True) + int(idx) -1
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

        """
        dup = self.copy()
        dup._label_flat()
        dup.data.sort_values(labels,axis=0,inplace=True,ascending=ascending)
        dup._label_condensed()
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

        """
        if set(self.shapes) & set(other.shapes) or self.data.shape[0] != other.data.shape[0]:
            raise KeyError('Dublicated keys or row count mismatch')
        else:
            dup = self.copy()
            dup.data = dup.data.join(other.data)
            for key in other.shapes:
                dup.shapes[key] = other.shapes[key]
            return dup


    def to_file(self,fname,format='ASCII',new_style=False):
        """
        Store as plain text file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file for writing.
        format : {ASCII'}, optional
            File format, defaults to 'ASCII'. Available formats are:
            - ASCII: Plain text file, extension '.txt'.
        new_style : Boolean, optional
            Write table in new style, indicating header lines by comment sign ('#') only.

        """
        def _to_ASCII(table,fname,new_style=False):
            """
            Store as plain text file.

            Parameters
            ----------
            table : Table object
                Table to write.
            fname : file, str, or pathlib.Path
                Filename or file for writing.
            new_style : Boolean, optional
                Write table in new style, indicating header lines by comment sign ('#') only.

            """
            seen = set()
            labels = []
            for l in [x for x in table.data.columns if not (x in seen or seen.add(x))]:
                if table.shapes[l] == (1,):
                    labels.append(f'{l}')
                elif len(table.shapes[l]) == 1:
                    labels += [f'{i+1}_{l}' \
                              for i in range(table.shapes[l][0])]
                else:
                    labels += [f'{util.srepr(table.shapes[l],"x")}:{i+1}_{l}' \
                              for i in range(np.prod(table.shapes[l]))]

            header = [f'# {comment}' for comment in table.comments] if new_style else \
                     [f'{len(table.comments)+1} header'] + table.comments

            try:
                f = open(fname,'w')
            except TypeError:
                f = fname

            for line in header + [' '.join(labels)]: f.write(line+'\n')
            table.data.to_csv(f,sep=' ',na_rep='nan',index=False,header=False)

        if format.lower() == 'ascii':
            return _to_ASCII(self,fname,new_style)
        else:
            raise TypeError(f'Unknown format {format}.')
