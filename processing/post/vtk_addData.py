#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,glob,re
import damask
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# -----------------------------
def findTag(filename,tag):
# -----------------------------
  
  with open(filename,'r') as myfile: 
    mypattern = re.compile(str(tag))
    for line in myfile:
      if mypattern.search(line): return True
  return False



# -----------------------------
# MAIN FUNCTION STARTS HERE
# -----------------------------

# --- input parsing

parser = OptionParser(usage='%prog [options] directory', description = """
Add data from an ASCII table to a VTK geometry file. 

""", version = scriptID)

parser.add_option('-s','--sub', action='store_true', dest='subdir', \
                  help='include files in subdirectories [%default]')
parser.set_defaults(subdir = False)

(options, dirname) = parser.parse_args()


# --- sanity checks

if dirname == []:
  parser.print_help()
  parser.error('no directory specified...')
else: 
  dirname = os.path.abspath(dirname[0])    # only use first argument

if not os.path.isdir(dirname):
  parser.print_help()
  parser.error('invalid directory "%s" specified...'%dirname)


# --- loop over "nodebased" and "ipbased" data files and 
#     copy data to corresponding geometry files

dataSetTag = {'nodebased':'POINT_DATA', 'ipbased':'CELL_DATA'}
for geomtype in ['nodebased','ipbased']:
  for vtkfilename in glob.iglob(dirname+os.sep+'*'+geomtype+'*.vtk'):

    if not os.path.dirname(vtkfilename) == dirname and not options.subdir: continue    # include files in subdir?
    datafilename = os.path.splitext(vtkfilename)[0] + '.txt'
    if not os.path.exists(datafilename): continue                                      # no corresponding datafile found 
    
    # --- read data from datafile

    with open(datafilename,'r') as datafile:                                           # open datafile in read mode
      table = damask.ASCIItable(fileIn=datafile)                                       # use ASCIItable class to read data file
      table.head_read()                                                                # read ASCII header info
      myData = []
      while table.data_read():                                                         # read line in datafile
        myData.append(table.data)
      myData = zip(*myData)                                                            # reorder data: first index now label, not node
    
    # --- append data to vtkfile 

    with open(vtkfilename,'a') as vtkfile:                                             # open vtkfile in append mode
      print vtkfilename
      if not findTag(vtkfilename,dataSetTag[geomtype]):                                # check if data set is already present...
        vtkfile.write(dataSetTag[geomtype] + ' %i'%len(myData[0]))                     # ... if not, write keyword
      for idx,label in enumerate(table.labels):                                        # write data
        vtkfile.write('\nSCALARS '+label+' double 1\nLOOKUP_TABLE default\n')          # all scalar data
        vtkfile.write('\n'.join(map(str,myData[idx])))
