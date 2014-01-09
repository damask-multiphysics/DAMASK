# $Id$

import numpy as np

try:
  import h5py
except:
  print('Could not import h5py.') 

class Result():
  '''
     General class for result parsing.
     Needs h5py to be installed
  '''
  
  def __init__(self,resultsFile):
    outFile=h5py.File(resultsFile,"a")
    print("Opened "+resultsFile+" with %i points"%outFile.attrs['Number of Materialpoints'])
    



