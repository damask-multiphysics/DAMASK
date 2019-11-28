import re

import pandas as pd
import numpy as np

class Table():
    """Store spreadsheet-like data."""
  
    def __init__(self,array,columns,comments=None):
        """
        New spreadsheet data.
        
        Parameters
        ----------
        array : numpy.ndarray
            Data.
        columns : dict
            Column labels and shape. Example 'F':(3,3) for a deformation gradient. 
        comments : iterable of str, optional
            Additional, human-readable information
        
        """
        self.data = pd.DataFrame(data=array)

        d = {}
        i = 0
        for label in columns:
            for components in range(np.prod(columns[label])):
                d[i] = label
                i+=1

        self.data.rename(columns=d,inplace=True)

        if comments is None:
            self.comments = []
        else:
            self.comments = [c for c in comments]

        self.columns = columns

    @staticmethod
    def from_ASCII(fname):
        """
        Create table from ASCII file.

        The first line needs to indicate the number of subsequent header lines as 'n header'.
        Vector data labels are indicated by '1_x, 2_x, ..., n_x'.
        Tensor data labels are indicated by '3x3:1_x, 3x3:2_x, ..., 3x3:9_x'.

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
       
        columns = {}
        for label in labels:
            tensor_column = re.search(r'[0-9,x]*?:[0-9]*?_',label)
            if tensor_column:
                my_shape = tensor_column.group().split(':',1)[0].split('x')
                columns[label.split('_',1)[1]] = tuple([int(d) for d in my_shape])
            else:
                vector_column = re.match(r'[0-9]*?_',label)
                if vector_column:
                    columns[label.split('_',1)[1]] = (int(label.split('_',1)[0]),)
                else:
                    columns[label]=(1,)
        
        return Table(np.loadtxt(f),columns,comments)

    def get_array(self,label):
        """Return data as array."""
        if re.match(r'[0-9]*?_',label):
            idx,key = label.split('_',1)
            return self.data[key].to_numpy()[:,int(idx)-1]
        else: 
            return self.data[label].to_numpy().reshape((-1,)+self.columns[label])

    def set_array(self,label,array,info):
        """
        Modify data in the spreadsheet.

        Parameters
        ----------
        label : str
            Label for the new data
        array : np.ndarray
            New data
        info : str
            Human-readable information about the new data

        """
        if np.prod(array.shape[1:],dtype=int) == 1: 
            self.comments.append('{}: {}'.format(label,info))
        else:
            self.comments.append('{} {}: {}'.format(label,array.shape[1:],info))

        if re.match(r'[0-9]*?_',label):
            idx,key = label.split('_',1)
            iloc = self.data.columns.get_loc(key).tolist().index(True) + int(idx) -1
            self.data.iloc[:,iloc] = array
        else: 
            self.data[label]       = array.reshape(self.data[label].shape)

    def get_labels(self):
        """Return the labels of all columns."""
        return [label for label in self.columns]

    def add_array(self,label,array,info):
        """
        Add data to the spreadsheet.
        
        Parameters
        ----------
        label : str
            Label for the new data
        array : np.ndarray
            New data
        info : str
            Human-readable information about the new data

        """
        if np.prod(array.shape[1:],dtype=int) == 1: 
            self.comments.append('{}: {}'.format(label,info))
        else:
            self.comments.append('{} {}: {}'.format(label,array.shape[1:],info))
        
        self.columns[label] = array.shape[1:] if len(array.shape) > 1 else (1,)
        size = np.prod(array.shape[1:],dtype=int) 
        new_data = pd.DataFrame(data=array.reshape(-1,size),
                                columns=[label for l in range(size)])
        self.data = pd.concat([self.data,new_data],axis=1)
        
    def to_ASCII(self,fname):
        """
        Store as plain text file.
        
        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file for reading.

        """

        labels = []
        for l in self.columns:
            if(self.columns[l] == (1,)):
                labels.append('{}'.format(l))
            elif(len(self.columns[l]) == 1):
                labels+=['{}_{}'.format(i+1,l)\
                         for i in range(self.columns[l][0])]
            else:
                labels+=['{}:{}_{}'.format(i+1,'x'.join([str(d) for d in self.columns[l]]),l)\
                               for i in range(np.prod(self.columns[l],dtype=int))]

        header = ['{} header'.format(len(self.comments)+1)]\
               + self.comments\
               + [' '.join(labels)]

        try:
            f = open(fname,'w')
        except TypeError:
            f = fname
        for line in header: f.write(line+'\n')
        self.data.to_csv(f,sep=' ',index=False,header=False)
