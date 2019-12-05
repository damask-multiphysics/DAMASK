import random
import re

import pandas as pd
import numpy as np

class Table():
    """Store spreadsheet-like data."""
  
    def __init__(self,data,shapes,comments=None):
        """
        New spreadsheet.
        
        Parameters
        ----------
        data : numpy.ndarray
            Data.
        shapes : dict with str:tuple pairs
            Shapes of the columns. Example 'F':(3,3) for a deformation gradient. 
        comments : iterable of str, optional
            Additional, human-readable information.
        
        """
        self.data = pd.DataFrame(data=data)

        labels = {}
        i = 0
        for label in shapes.keys():
            for components in range(np.prod(shapes[label])):
                labels[i] = label
                i+=1

        if i != self.data.shape[1]:
            raise IndexError('Shape mismatch between shapes and data')

        self.data.rename(columns=labels,inplace=True)

        if comments is None:
            self.comments = []
        else:
            self.comments = [c for c in comments]

        self.shapes = shapes

    @staticmethod
    def from_ASCII(fname):
        """
        Create table from ASCII file.

        The first line needs to indicate the number of subsequent header lines as 'n header'.
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

        header,keyword = f.readline().split()
        if keyword == 'header':
            header = int(header)
        else:
            raise Exception
        comments = [f.readline()[:-1] for i in range(header-1)]
        labels   = f.readline().split()
       
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
                    shapes[label]=(1,)
        
        data = pd.read_csv(f,names=[i for i in range(len(labels))],sep=r'\s+').to_numpy()
        
        return Table(data,shapes,comments)


    def labels(self):
        """Return the labels of all columns."""
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
            return self.data[key].to_numpy()[:,int(idx)-1].reshape((-1,1))
        else: 
            return self.data[label].to_numpy().reshape((-1,)+self.shapes[label])

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
        if info is not None:
            if np.prod(data.shape[1:],dtype=int) == 1: 
                self.comments.append('{}: {}'.format(label,info))
            else:
                self.comments.append('{} {}: {}'.format(label,data.shape[1:],info))

        if re.match(r'[0-9]*?_',label):
            idx,key = label.split('_',1)
            iloc = self.data.columns.get_loc(key).tolist().index(True) + int(idx) -1
            self.data.iloc[:,iloc] = data
        else: 
            self.data[label]       = data.reshape(self.data[label].shape)

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
        if info is not None:
           if np.prod(data.shape[1:],dtype=int) == 1: 
               self.comments.append('{}: {}'.format(label,info))
           else:
               self.comments.append('{} {}: {}'.format(label,data.shape[1:],info))
        
        self.shapes[label] = data.shape[1:] if len(data.shape) > 1 else (1,)
        size = np.prod(data.shape[1:],dtype=int) 
        new_data = pd.DataFrame(data=data.reshape(-1,size),
                                columns=[label for l in range(size)])
        self.data = pd.concat([self.data,new_data],axis=1)

    def delete(self,label):
        """
        Delete column data.

        Parameters
        ----------
        label : str
            Column label.

        """
        self.data.drop(columns=label,inplace=True)

        del self.shapes[label]

    def rename(self,label_old,label_new,info=None):
        """
        Rename column data.

        Parameters
        ----------
        label_old : str
            Old column label.
        label_new : str
            New column label.

        """
        self.data.rename(columns={label_old:label_new},inplace=True)

        comments = '{} => {}'.format(label_old,label_new)
        comments += ': {}'.format(info) if info is not None else ''
        self.comments.append(comments)

        self.shapes[label_new] = self.shapes.pop(label_old)


    def sort_by(self,labels,ascending=True):
        """
        Get column data.

        Parameters
        ----------
        label : list of str or str
            Column labels.
        ascending : bool, optional
            Set sort order.

        """
        _temp = []
        _labels = []
        for label in labels if isinstance(labels,list) else [labels]:
            if re.match(r'[0-9]*?_',label):
                _temp.append(str(random.getrandbits(128)))
                self.add(_temp[-1],self.get(label))
                _labels.append(_temp[-1])
            else: 
                _labels.append(label)

        self.data.sort_values(_labels,axis=0,inplace=True,ascending=ascending)
        for t in _temp: self.delete(t)
        self.comments.append('sorted by [{}]'.format(', '.join(labels)))

    def to_ASCII(self,fname):
        """
        Store as plain text file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file for reading.

        """
        labels = []
        for l in self.shapes:
            if(self.shapes[l] == (1,)):
                labels.append('{}'.format(l))
            elif(len(self.shapes[l]) == 1):
                labels+=['{}_{}'.format(i+1,l)\
                         for i in range(self.shapes[l][0])]
            else:
                labels+=['{}:{}_{}'.format('x'.join([str(d) for d in self.shapes[l]]),i+1,l)\
                               for i in range(np.prod(self.shapes[l],dtype=int))]

        header = ['{} header'.format(len(self.comments)+1)]\
               + self.comments\
               + [' '.join(labels)]

        try:
            f = open(fname,'w')
        except TypeError:
            f = fname
        for line in header: f.write(line+'\n')
        self.data.to_csv(f,sep=' ',index=False,header=False)
