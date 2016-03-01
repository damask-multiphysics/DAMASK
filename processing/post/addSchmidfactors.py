#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

slipnormal_temp = [
    [0,0,0,1],
    [0,0,0,1],
    [0,0,0,1],
    [0,1,-1,0],
    [-1,0,1,0],
    [1,-1,0,0],
    [0,1,-1,1],
    [-1,1,0,1],
    [-1,0,1,1],
    [0,-1,1,1],
    [1,-1,0,1],
    [1,0,-1,1],
    [0,1,-1,1],
    [0,1,-1,1],
    [-1,1,0,1],
    [-1,1,0,1],
    [-1,0,1,1],
    [-1,0,1,1],
    [0,-1,1,1],
    [0,-1,1,1],
    [1,-1,0,1],
    [1,-1,0,1],
    [1,0,-1,1],
    [1,0,-1,1],
   ]

slipdirection_temp = [
    [2,-1,-1,0],
    [-1,2,-1,0],
    [-1,-1,2,0],
    [2,-1,-1,0],
    [-1,2,-1,0],
    [-1,-1,2,0],
    [2,-1,-1,0],
    [1,1,-2,0],
    [-1,2,-1,0],
    [-2,1,1,0],
    [-1,-1,2,0],
    [1,-2,1,0],
    [-1,2,-1,3],
    [1,1,-2,3],
    [-2,1,1,3],
    [-1,2,-1,3],
    [-1,-1,2,3],
    [-2,1,1,3],
    [1,-2,1,3],
    [-1,-1,2,3],
    [2,-1,-1,3],
    [1,-2,1,3],
    [1,1,-2,3],
    [2,-1,-1,3],
   ]

# slip normals and directions according to cpfem implementation
Nslipsystems = {'fcc': 12, 'bcc': 24, 'hex': 24}
slipnormal = { \
   'fcc': [
    [1,1,1],
    [1,1,1],
    [1,1,1],
    [-1,-1,1],
    [-1,-1,1],
    [-1,-1,1],
    [1,-1,-1],
    [1,-1,-1],
    [1,-1,-1],
    [-1,1,-1],
    [-1,1,-1],
    [-1,1,-1],
   ],
   'bcc': [
    [0,1,1],
    [0,1,1],
    [0,-1,1],
    [0,-1,1],
    [1,0,1],
    [1,0,1],
    [-1,0,1],
    [-1,0,1],
    [1,1,0],
    [1,1,0],
    [-1,1,0],
    [-1,1,0],
    [2,1,1],
    [-2,1,1],
    [2,-1,1],
    [2,1,-1],
    [1,2,1],
    [-1,2,1],
    [1,-2,1],
    [1,2,-1],
    [1,1,2],
    [-1,1,2],
    [1,-1,2],
    [1,1,-2],
   ],
   'hex': [  # these are dummy numbers and are recalculated based on the above hex real slip systems.
    [1,1,0],
    [1,1,0],
    [1,0,1],
    [1,0,1],
    [0,1,1],
    [0,1,1],
    [1,-1,0],
    [1,-1,0],
    [-1,0,1],
    [-1,0,1],
    [0,-1,1],
    [0,-1,1],
    [2,-1,1],
    [1,-2,-1],
    [1,1,2],
    [2,1,1],
    [1,2,-1],
    [1,-1,2],
    [2,1,-1],
    [1,2,1],
    [1,-1,-2],
    [2,-1,-1],
    [1,-2,1],
    [1,1,-2],
   ],   
   }
slipdirection = { \
   'fcc': [
    [0,1,-1],
    [-1,0,1],
    [1,-1,0],
    [0,-1,-1],
    [1,0,1],
    [-1,1,0],
    [0,-1,1],
    [-1,0,-1],
    [1,1,0],
    [0,1,1],
    [1,0,-1],
    [-1,-1,0],
    ],
   'bcc': [
    [1,-1,1],
    [-1,-1,1],
    [1,1,1],
    [-1,1,1],
    [-1,1,1],
    [-1,-1,1],
    [1,1,1],
    [1,-1,1],
    [-1,1,1],
    [-1,1,-1],
    [1,1,1],
    [1,1,-1],
    [-1,1,1],
    [1,1,1],
    [1,1,-1],
    [1,-1,1],
    [1,-1,1],
    [1,1,-1],
    [1,1,1],
    [-1,1,1],
    [1,1,-1],
    [1,-1,1],
    [-1,1,1],
    [1,1,1],
   ],
   'hex': [ # these are dummy numbers and are recalculated based on the above hex real slip systems.
    [-1,1,1],
    [1,-1,1],
    [-1,-1,1],
    [-1,1,1],
    [-1,-1,1],
    [1,-1,1],
    [1,1,1],
    [-1,-1,1],
    [1,-1,1],
    [1,1,1],
    [1,1,1],
    [-1,1,1],
    [1,1,-1],
    [1,1,-1],
    [1,1,-1],
    [1,-1,-1],
    [1,-1,-1],
    [1,-1,-1],
    [1,-1,1],
    [1,-1,1],
    [1,-1,1],
    [1,1,1],
    [1,1,1],
    [1,1,1],
   ],
   }

def applyEulers(phi1,Phi,phi2,x):
  """transform x given in crystal coordinates to xbar returned in lab coordinates for Euler angles phi1,Phi,phi2"""
  eulerRot = [[ math.cos(phi1)*math.cos(phi2) - math.cos(Phi)*math.sin(phi1)*math.sin(phi2),
               -math.cos(phi1)*math.sin(phi2) - math.cos(Phi)*math.cos(phi2)*math.sin(phi1),
                math.sin(Phi)*math.sin(phi1)
              ],
              [ math.cos(phi2)*math.sin(phi1) + math.cos(Phi)*math.cos(phi1)*math.sin(phi2),
                math.cos(Phi)*math.cos(phi1)*math.cos(phi2) - math.sin(phi1)*math.sin(phi2),
               -math.sin(Phi)*math.cos(phi1)
              ],
              [ math.sin(Phi)*math.sin(phi2),
                math.sin(Phi)*math.cos(phi2),
                math.cos(Phi)
             ]]
  
  xbar = [0,0,0]
  if len(x) == 3:
    for i in range(3): 
      xbar[i] = sum([eulerRot[i][j]*x[j] for j in range(3)])
  return xbar

def normalize(x):
  
  norm = math.sqrt(sum([x[i]*x[i] for i in range(len(x))]))
  
  return [x[i]/norm for i in range(len(x))]

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add columns listing Schmid factors (and optional trace vector of selected system) for given Euler angles.

""", version = scriptID)

parser.add_option('-l','--lattice',
                  dest='lattice', type='choice', choices=('fcc','bcc','hex'), metavar='string',
                  help="type of lattice structure [%default] {fcc,bcc',hex}")
parser.add_option('--direction',
                  dest='forcedirection', type='int', nargs=3, metavar='int int int',
                  help='force direction in lab coordinates %default')
parser.add_option('-n','--normal',
                  dest='stressnormal', type='int', nargs=3, metavar='int int int',
                  help='stress plane normal in lab coordinates ')
parser.add_option('--trace',
                  dest='traceplane', type='int', nargs=3, metavar='int int int',
                  help='normal (in lab coordinates) of plane on which the plane trace of the Schmid factor(s) is reported')
parser.add_option('--covera',
                  dest='CoverA', type='float', metavar='float',
                  help='C over A ratio for hexagonal systems')
parser.add_option('-r','--rank',
                  dest='rank', type='int', nargs=3, metavar='int int int',
                  help="report trace of r'th highest Schmid factor [%default]")
parser.add_option('-e', '--eulers',
                  dest='eulers', metavar='string',
                  help='Euler angles label')
parser.add_option('-d', '--degrees',
                  dest='degrees', action='store_true',
                  help='Euler angles are given in degrees [%default]')
parser.set_defaults(lattice = 'fcc')
parser.set_defaults(forcedirection = [0, 0, 1])
parser.set_defaults(stressnormal = None)
parser.set_defaults(traceplane = None)
parser.set_defaults(rank = 0)
parser.set_defaults(CoverA = 1.587)
parser.set_defaults(eulers = 'eulerangles')

(options,filenames) = parser.parse_args()

options.forcedirection = normalize(options.forcedirection)
if options.stressnormal:
  if abs(sum([options.forcedirection[i] * options.stressnormal[i] for i in range(3)])) < 1e-3:
    options.stressnormal = normalize(options.stressnormal)
  else:
    parser.error('stress plane normal not orthogonal to force direction')
else:
  options.stressnormal = options.forcedirection
if options.traceplane:
  options.traceplane = normalize(options.traceplane)
options.rank = min(options.rank,Nslipsystems[options.lattice])

datainfo = {                                                                                        # list of requested labels per datatype
             'vector':     {'len':3,
                            'label':[]},
           }

datainfo['vector']['label'] += [options.eulers]

toRadians = math.pi/180.0 if options.degrees else 1.0                                               # rescale degrees to radians
# Convert 4 Miller indices notation of hex to orthogonal 3 Miller indices notation
if options.lattice=='hex':                                                                          
  for i in range(Nslipsystems[options.lattice]):
    slipnormal[options.lattice][i][0]=slipnormal_temp[i][0]
    slipnormal[options.lattice][i][1]=(slipnormal_temp[i][0]+2.0*slipnormal_temp[i][1])/math.sqrt(3.0)
    slipnormal[options.lattice][i][2]=slipnormal_temp[i][3]/options.CoverA
    slipdirection[options.lattice][i][0]=slipdirection_temp[i][0]*1.5                               # direction [uvtw]->[3u/2 (u+2v)*sqrt(3)/2 w*(c/a)] ,
    slipdirection[options.lattice][i][1]=(slipdirection_temp[i][0]+2.0*slipdirection_temp[i][1])*(0.5*math.sqrt(3.0))
    slipdirection[options.lattice][i][2]=slipdirection_temp[i][3]*options.CoverA

  for i in range(Nslipsystems[options.lattice]):
    slipnormal[options.lattice][i]=normalize(slipnormal[options.lattice][i])
    slipdirection[options.lattice][i]=normalize(slipdirection[options.lattice][i])

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,buffered = False)
  except:
    continue
  damask.util.report(scriptName,name)

  table.head_read()                                                                                 # read ASCII header info
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))

  key = '1_%s'%datainfo['vector']['label'][0]
  if key not in table.labels:
    file['croak'].write('column %s not found...\n'%key)
    continue
  else:
    column = table.labels.index(key)                                                                # remember columns of requested data

# ------------------------------------------ assemble header ---------------------------------------

  table.labels_append(['%i_S(%i_%i_%i)[%i_%i_%i]'%(i+1,
                       slipnormal[options.lattice][i][0],
                       slipnormal[options.lattice][i][1],
                       slipnormal[options.lattice][i][2],
                       slipdirection[options.lattice][i][0],
                       slipdirection[options.lattice][i][1],
                       slipdirection[options.lattice][i][2],
                     ) for i in range(Nslipsystems[options.lattice])])

  if options.traceplane:
    if options.rank > 0:
      table.labels_append('trace_x trace_y trace_z system')
    else:
      table.labels_append(['(%i)tx\tty\ttz'%(i+1) for i in range(Nslipsystems[options.lattice])])
  table.head_write()

# ------------------------------------------ process data ------------------------------------------
  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    [phi1,Phi,phi2] = Eulers=toRadians*np.array(map(\
                           float,table.data[column:column+datainfo['vector']['len']]))
    S = [ sum( [applyEulers(phi1,Phi,phi2,normalize( \
                slipnormal[options.lattice][slipsystem]))[i]*options.stressnormal[i] for i in range(3)] ) * \
          sum( [applyEulers(phi1,Phi,phi2,normalize( \
                slipdirection[options.lattice][slipsystem]))[i]*options.forcedirection[i] for i in range(3)] ) \
          for slipsystem in range(Nslipsystems[options.lattice]) ]
    table.data_append(S)
    if options.traceplane:
      trace = [np.cross(options.traceplane,applyEulers(phi1,Phi,phi2,normalize(slipnormal[options.lattice][slipsystem]))) \
                 for slipsystem in range(Nslipsystems[options.lattice]) ]
      if options.rank == 0:
        table.data_append('\t'.join(map(lambda x:'%f\t%f\t%f'%(x[0],x[1],x[2]),trace)))
      elif options.rank > 0:
        SabsSorted = sorted([(abs(S[i]),i) for i in range(len(S))])
        table.data_append('\t'.join(map(str,trace[SabsSorted[-options.rank][1]])) + '\t%i'%(1+SabsSorted[-options.rank][1]))
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)
