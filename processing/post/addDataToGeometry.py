#!/usr/bin/env python 

import os, sys, string, glob
import damask
from optparse import OptionParser



# -----------------------------
def writeHeader(myfile,stat,geomtype):
# -----------------------------
    
  myfile.write('2\theader\n')
  myfile.write(string.replace('$Id: $','\n','\\n')+
           '\t' + ' '.join(sys.argv[1:]) + '\n')
  if geomtype == 'nodebased':
    myfile.write('node')
    for i in range(stat['NumberOfNodalScalars']):
      myfile.write('\t%s'%''.join(stat['LabelOfNodalScalar'][i].split()))
    
  elif geomtype == 'ipbased':
    myfile.write('elem\tip')
    for i in range(stat['NumberOfElementalScalars']):
      myfile.write('\t%s'%''.join(stat['LabelOfElementalScalar'][i].split()))
  
  myfile.write('\n')
   
  return True



# -----------------------------
# MAIN FUNCTION STARTS HERE
# -----------------------------

# --- input parsing

parser = OptionParser(usage='%prog [options] directory', description = """
Add data from an ASCII table to a VTK geometry file. 
""" + string.replace('$Id:  $','\n','\\n')
)

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

for geomtype in ['nodebased','ipbased']:
  for vtkfilename in glob.iglob(dirname+os.sep+'*'+geomtype+'*.vtk'):
    
    if not os.path.dirname(vtkfilename) == dirname and not options.subdir: continue  # include files in subdir?
    datafilename = os.path.splitext(vtkfilename)[0] + '.txt'
    if not os.path.exists(datafilename): continue                                    # no corresponding datafile found 
    
    with open(vtkfilename,'a') as vtkfile:
      print vtkfilename
      with open(datafilename,'r') as datafile:

        table = damask.ASCIItable(fileIn=datafile)                                   # use ASCIItable class to read data file
        table.head_read()                                                            # read ASCII header info
        myData = []
        while table.data_read():                                                     # read line in datafile
          myData.append(table.data)
        myData = zip(*myData)                                                        # reorder data: first index now label, not node
        vtkfile.write('CELL_DATA %i'%len(myData[0]))
        for idx,label in enumerate(table.labels):
          vtkfile.write('\nSCALARS '+label+' float 1\nLOOKUP_TABLE default\n')       # all scalar data
          vtkfile.write('\n'.join(map(str,myData[idx])))


# ---------------------------       DONE     --------------------------------
