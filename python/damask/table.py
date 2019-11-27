import re

import pandas as pd
import numpy as np

class Table():
    """Store spreadsheet-like data."""
  
    def __init__(self,array,headings,comments=None):
        """
        New spreadsheet data.
        
        Parameters
        ----------
        array : numpy.ndarray
            Data.
        headings : dict
            Column headings. Labels as keys and shape as tuple. Example 'F':(3,3) for a deformation gradient. 
        comments : iterable of str, optional
            Additional, human-readable information
        
        """
        self.data = pd.DataFrame(data=array)

        d = {}
        i = 0
        for label in headings:
            for components in range(np.prod(headings[label])):
                d[i] = label
                i+=1

        self.data.rename(columns=d,inplace=True)

        if comments is None:
            self.comments = []
        else:
            self.comments = [c for c in comments]

        self.headings  = headings

    @staticmethod
    def from_ASCII(fname):
        """
        Create table from ASCII file.

        The first line needs to indicate the number of subsequent header lines as 'n header'.
        Vector data labels are indicated by '1_x, 2_x, ..., n_x'.
        Tensor data labels are indicated by '3x3:1_x, 3x3:2_x, ..., 3x3:9_x'.
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
       
        headings = {}
        for label in labels:
            tensor_column = re.search(r'[0-9,x]*?:[0-9]*?_',label)
            if tensor_column:
                my_shape = tensor_column.group().split(':',1)[0].split('x')
                headings[label.split('_',1)[1]] = tuple([int(d) for d in my_shape])
            else:
                vector_column = re.match(r'[0-9]*?_',label)
                if vector_column:
                    headings[label.split('_',1)[1]] = (int(label.split('_',1)[0]),)
                else:
                    headings[label]=(1,)
        
        return Table(np.loadtxt(f),headings,comments)

    def get_array(self,label):
        """Return data as array."""
        if re.match(r'[0-9]*?_',label):
            idx,key = label.split('_',1)
            return self.data[key].to_numpy()[:,int(idx)-1]
        else: 
            return self.data[label].to_numpy().reshape((-1,)+self.headings[label])

    def set_array(self,label,array):
        """Set data."""
        if re.match(r'[0-9]*?_',label):
            idx,key = label.split('_',1)
            iloc = self.data.columns.get_loc(key).tolist().index(True) + int(idx) -1
            self.data.iloc[:,iloc] = array
        else: 
            self.data[label]       = array


    def get_labels(self):
        """Return the labels of all columns."""
        return [label for label in self.headings]

    def add_array(self,label,array,info):
        if np.prod(array.shape[1:],dtype=int) == 1: 
            self.comments.append('{}: {}'.format(label,info))
        else:
            self.comments.append('{} {}: {}'.format(label,array.shape[1:],info))
        
        self.headings[label] = array.shape[1:] if len(array.shape) > 1 else (1,)
        size = np.prod(array.shape[1:],dtype=int) 
        new_data = pd.DataFrame(data=array.reshape(-1,size),
                                columns=[label for l in range(size)])
        self.data = pd.concat([self.data,new_data],axis=1)
        
    def to_ASCII(self,fname):
        labels = []
        for l in self.headings:
            if(self.headings[l] == (1,)):
                labels.append('{}'.format(l))
            elif(len(self.headings[l]) == 1):
                labels+=['{}_{}'.format(i+1,l)\
                         for i in range(self.headings[l][0])]
            else:
                labels+=['{}:{}_{}'.format(i+1,'x'.join([str(d) for d in self.headings[l]]),l)\
                               for i in range(np.prod(self.headings[l],dtype=int))]

        header = ['{} header'.format(len(self.comments)+1)]\
               + self.comments\
               + [' '.join(labels)]

        try:
            f = open(fname,'w')
        except TypeError:
            f = fname
        for line in header: f.write(line+'\n')
        self.data.to_csv(f,sep=' ',index=False,header=False)
