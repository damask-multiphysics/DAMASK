import re

import pandas as pd
import numpy as np

class Table():
    """Read and write to ASCII tables"""
    
    def __init__(self,name):
        self.name = name
        with open(self.name) as f:
            header,keyword = f.readline().split()
            if keyword == 'header':
                header = int(header)
            else:
                raise Exception
            self.comments = [f.readline()[:-1] for i in range(header-1)]
            labels_raw    = f.readline().split()
            self.data     = pd.read_csv(f,delim_whitespace=True,header=None)
            
            labels_repeated = [l.split('_',1)[1] if '_' in l else l for l in labels_raw]
            self.data.rename(columns=dict(zip([l for l in self.data.columns],labels_repeated)),inplace=True)
            
            self.shape = {}
            for l in labels_raw:
                tensor_column = re.search(':.*?_',l)
                if tensor_column:
                    my_shape = tensor_column.group()[1:-1].split('x')
                    self.shape[l.split('_',1)[1]] = tuple([int(d) for d in my_shape])
                else:
                    vector_column = re.match('.*?_',l)
                    if vector_column:
                        self.shape[l.split('_',1)[1]] = (int(l.split('_',1)[0]),)
                    else:
                        self.shape[l]=(1,)
            
            self.labels = list(dict.fromkeys(labels_repeated))


    def get_array(self,label):
        return self.data[label].to_numpy().reshape((-1,)+self.shape[label])


    def add_array(self,label,array,info):
        if np.product(array.shape[1:],dtype=int) == 1: 
            self.comments.append('{}: {}'.format(label,info))

        else:
            self.comments.append('{} {}: {}'.format(label,array.shape[1:],info))
        
        self.shape[label] = array.shape[1:]
        self.labels.append(label)
        size = np.product(array.shape[1:]) 
        new_data = pd.DataFrame(data=array.reshape(-1,size),
                                columns=[label for l in range(size)])
        self.data = pd.concat([self.data,new_data],axis=1)
        

    def to_ASCII(self,name=None):
        labels = []
        for l in self.labels:
            if(self.shape[l] == (1,)):
                labels.append('{}'.format(l))
            elif(len(self.shape[l]) == 1):
                labels+=['{}_{}'.format(i+1,l)\
                         for i in range(self.shape[l][0])]
            else:
                labels+=['{}:{}_{}'.format(i+1,'x'.join([str(d) for d in self.shape[l]]),l)\
                               for i in range(np.product(self.shape[l]))]

        header = ['{} header'.format(len(self.comments)+1)]\
               + self.comments\
               + [' '.join(labels)]

        with open(name if name is not None else self.name,'w') as f:
            for line in header: f.write(line+'\n')
            self.data.to_csv(f,sep=' ',index=False,header=False)
