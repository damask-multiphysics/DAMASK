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
        try:
            f = open(fname)
        except TypeError:
            f = fname

        header,keyword = f.readline().split()
        if keyword == 'header':
            header = int(header)
        else:
            raise Exception
        comments   = [f.readline()[:-1] for i in range(header-1)]
        labels_raw = f.readline().split()
        labels     = [l.split('_',1)[1] if '_' in l else l for l in labels_raw]
       
        headings = {}
        for l in labels_raw:
            tensor_column = re.search(':.*?_',l)
            if tensor_column:
                my_shape = tensor_column.group()[1:-1].split('x')
                headings[l.split('_',1)[1]] = tuple([int(d) for d in my_shape])
            else:
                vector_column = re.match('.*?_',l)
                if vector_column:
                    headings[l.split('_',1)[1]] = (int(l.split('_',1)[0]),)
                else:
                    headings[l]=(1,)
        
        return Table(np.loadtxt(f),headings,comments)


    def get_array(self,label):
        return self.data[label].to_numpy().reshape((-1,)+self.headings[label])


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
