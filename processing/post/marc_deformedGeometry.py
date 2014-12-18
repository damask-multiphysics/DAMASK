#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os, sys, math, string, numpy, shutil
import damask
from optparse import OptionParser

scriptID = '$Id$'
scriptName = os.path.splitext(scriptID.split()[1])[0]

# -----------------------------
# MAIN FUNCTION STARTS HERE
# -----------------------------

# --- input parsing

parser = OptionParser(usage='%prog [options] resultfile', description = """
Create vtk files for the (deformed) geometry that belongs to a .t16 (MSC.Marc) results file.

""", version = scriptID)

parser.add_option('-d','--dir', dest='dir', \
                  help='name of subdirectory to hold output [%default]')
parser.add_option('-r','--range', dest='range', type='int', nargs=3, \
                  help='range of positions (or increments) to output (start, end, step) [all]')
parser.add_option('--increments', action='store_true', dest='getIncrements', \
                  help='switch to increment range [%default]')
parser.add_option('-t','--type', dest='type', type='choice', choices=['ipbased','nodebased'], \
                  help='processed geometry type [ipbased and nodebased]')

parser.set_defaults(dir = 'vtk')
parser.set_defaults(getIncrements= False)

(options, files) = parser.parse_args()

# --- basic sanity checks

if files == []:
  parser.print_help()
  parser.error('no file specified...')

filename = os.path.splitext(files[0])[0]
if not os.path.exists(filename+'.t16'):
  parser.print_help()
  parser.error('invalid file "%s" specified...'%filename+'.t16')

if not options.type :
  options.type = ['nodebased', 'ipbased']
else: 
  options.type = [options.type]


# --- more sanity checks

sys.path.append(damask.solver.Marc().libraryPath('../../'))
try:
  from py_post import *
except:
  print('error: no valid Mentat release found')
  sys.exit(-1)


# ---------------------------   open results file and initialize mesh    ----------
    
p = post_open(filename+'.t16')
p.moveto(0)
Nnodes = p.nodes()
Nincrements = p.increments() - 1                   # t16 contains one "virtual" increment (at 0)
if damask.core.mesh.mesh_init_postprocessing(filename+'.mesh') > 0:
  print('error: init not successful')
  sys.exit(-1)
Ncellnodes = damask.core.mesh.mesh_get_Ncellnodes()
unitlength = damask.core.mesh.mesh_get_unitlength()


# ---------------------------   create output dir   --------------------------------

dirname = os.path.abspath(os.path.join(os.path.dirname(filename),options.dir))
if not os.path.isdir(dirname):
  os.mkdir(dirname,0755)


# ---------------------------   get positions   --------------------------------

incAtPosition = {}
positionOfInc = {}

for position in range(Nincrements):
  p.moveto(position+1)
  incAtPosition[position] = p.increment            # remember "real" increment at this position
  positionOfInc[p.increment] = position            # remember position of "real" increment

if not options.range:
  options.getIncrements = False
  locations = range(Nincrements)                   # process all positions
else:
  options.range = list(options.range)              # convert to list
  if options.getIncrements:
    locations = [positionOfInc[x] for x in range(options.range[0],options.range[1]+1,options.range[2])
                                   if x in positionOfInc]
  else:
    locations = range( max(0,options.range[0]),
                       min(Nincrements,options.range[1]+1),
                       options.range[2] )

increments = [incAtPosition[x] for x in locations] # build list of increments to process



# ---------------------------   loop over positions   --------------------------------

for incCount,position in enumerate(locations):     # walk through locations

  p.moveto(position+1)                             # wind to correct position

  # --- get displacements 

  node_displacement = [[0,0,0] for i in range(Nnodes)]
  for n in range(Nnodes):
    if p.node_displacements():
      node_displacement[n] = map(lambda x:x*unitlength,list(p.node_displacement(n)))
  c = damask.core.mesh.mesh_build_cellnodes(numpy.array(node_displacement).T,Ncellnodes)
  cellnode_displacement = [[c[i][n] for i in range(3)] for n in range(Ncellnodes)]


  # --- append displacements to corresponding files
  
  for geomtype in options.type:
    outFilename = eval('"'+eval("'%%s_%%s_inc%%0%ii.vtk'%(math.log10(max(increments+[1]))+1)")+'"%(dirname + os.sep + os.path.split(filename)[1],geomtype,increments[incCount])')
    print outFilename
    shutil.copyfile('%s_%s.vtk'%(filename,geomtype),outFilename)
  
    with open(outFilename,'a') as myfile:
      myfile.write("POINT_DATA %i\n"%{'nodebased':Nnodes,'ipbased':Ncellnodes}[geomtype])
      myfile.write("VECTORS displacement double\n")
      coordinates = {'nodebased':node_displacement,'ipbased':cellnode_displacement}[geomtype]
      for n in range({'nodebased':Nnodes,'ipbased':Ncellnodes}[geomtype]):
        myfile.write("%.8e %.8e %.8e\n"%(coordinates[n][0],coordinates[n][1],coordinates[n][2]))



# ---------------------------       DONE     --------------------------------
