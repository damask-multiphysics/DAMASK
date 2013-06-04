#!/usr/bin/python

import os,re,sys,math
from optparse import OptionParser

CoverA=1.587
slipnormal_temp = [ # This is the real slip system information for hex aka titanium for now.
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


# --------------------------------------------------------------------

def applyEulers(phi1,Phi,phi2,x):
  """ transform x given in crystal coordinates to xbar returned in lab coordinates for Euler angles phi1,Phi,phi2 """
  
  eulerRot = [[ math.cos(phi1)*math.cos(phi2) - math.cos(Phi)*math.sin(phi1)*math.sin(phi2), - math.cos(phi1)*math.sin(phi2) - math.cos(Phi)*math.cos(phi2)*math.sin(phi1),  math.sin(Phi)*math.sin(phi1)], \
              [ math.cos(phi2)*math.sin(phi1) + math.cos(Phi)*math.cos(phi1)*math.sin(phi2),   math.cos(Phi)*math.cos(phi1)*math.cos(phi2) - math.sin(phi1)*math.sin(phi2), -math.sin(Phi)*math.cos(phi1)], \
              [                                 math.sin(Phi)*math.sin(phi2),                                   math.sin(Phi)*math.cos(phi2),            math.cos(Phi)]]
  
  xbar = [0,0,0]
  if len(x) == 3:
    for i in range(3): 
      xbar[i] = sum([eulerRot[i][j]*x[j] for j in range(3)])
  return xbar

# --------------------------------------------------------------------

def normalize(x):
  
  norm = math.sqrt(sum([x[i]*x[i] for i in range(len(x))]))
  
  return [x[i]/norm for i in range(len(x))]

# --------------------------------------------------------------------

def crossproduct(x,y):
  
  return [
           x[1]*y[2]-y[1]*x[2],
           x[2]*y[0]-y[2]*x[0],
           x[0]*y[1]-y[0]*x[1],
         ]

# --------------------------------------------------------------------



# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
parser = OptionParser(usage='%prog [options] [file]', description = """
Add columns listing Schmid factors (and optional trace vector of selected system) for given Euler angles.
Column headings need to have names 'phi1', 'Phi', 'phi2'.

$Id$
""")

parser.add_option('-l','--lattice', dest='lattice', choices=('fcc','bcc','hex'), \
                                    help='key for lattice type [%default]')
parser.add_option('-d','--forcedirection', dest='forcedirection', type='int', nargs=3, \
                                      help='force direction in lab coordinates [%default]')
parser.add_option('-n','--stressnormal', dest='stressnormal', type='int', nargs=3, \
                                      help='stress plane normal in lab coordinates [%default]')
parser.add_option('-t','--trace', dest='traceplane', type='int', nargs=3, \
                                      help="normal (in lab coordinates) of plane on which the plane trace of the Schmid factor(s) is reported [%default]")
parser.add_option('-r','--rank', dest='rank', type='int', \
                                      help="report trace of r'th highest Schmid factor [%default]")
parser.set_defaults(lattice = 'fcc')
parser.set_defaults(forcedirection = [0, 0, 1])
parser.set_defaults(stressnormal = None)
parser.set_defaults(traceplane = None)
parser.set_defaults(rank = 0)

(options,filename) = parser.parse_args()

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

# read from standard input unless input file specified
if filename == []:
  file = sys.stdin
elif os.path.exists(filename[0]):
  file = open(filename[0])

# read data
content = file.readlines()
file.close()

#  get labels by either read the first row, or - if keyword header is present - the last line of the header
headerlines = 1
m = re.search('(\d+)\s*head', content[0].lower())
if m:
  headerlines = int(m.group(1))+1
labels = content[headerlines-1].split()
data = content[headerlines:]

# Convert 4 Miller indices notation of hex to orthogonal 3 Miller indices notation
if options.lattice=="hex":
  for i in range(Nslipsystems[options.lattice]):
    slipnormal[options.lattice][i][0]=slipnormal_temp[i][0]
    slipnormal[options.lattice][i][1]=(slipnormal_temp[i][0]+2.0*slipnormal_temp[i][1])/math.sqrt(3.0)
    slipnormal[options.lattice][i][2]=slipnormal_temp[i][3]/CoverA
    slipdirection[options.lattice][i][0]=slipdirection_temp[i][0]*1.5 # direction [uvtw]->[3u/2 (u+2v)*sqrt(3)/2 w*(c/a)] ,
    slipdirection[options.lattice][i][1]=(slipdirection_temp[i][0]+2.0*slipdirection_temp[i][1])*(0.5*math.sqrt(3.0))
    slipdirection[options.lattice][i][2]=slipdirection_temp[i][3]*CoverA

  for i in range(Nslipsystems[options.lattice]):
    slipnormal[options.lattice][i]=normalize(slipnormal[options.lattice][i])
    slipdirection[options.lattice][i]=normalize(slipdirection[options.lattice][i])

for c in range(len(labels)):
  m = re.search('.*([Pp]hi\d*).*', labels[c])
  if m:
    if m.group(1).lower() == "phi1":
      phi1Column = c
    elif m.group(1).lower() == "phi":
      PhiColumn = c
    elif m.group(1).lower() == "phi2":
      phi2Column = c

output = '1\theader\n' + \
         '\t'.join(map(str,labels)) + \
         '\t' + \
         '\t'.join(['(%i)S(%i %i %i)[%i %i %i]'%(i+1,
                                             slipnormal[options.lattice][i][0],
                                             slipnormal[options.lattice][i][1],
                                             slipnormal[options.lattice][i][2],
                                             slipdirection[options.lattice][i][0],
                                             slipdirection[options.lattice][i][1],
                                             slipdirection[options.lattice][i][2],
                                            ) for i in range(Nslipsystems[options.lattice])])
if options.traceplane:
  if options.rank > 0:
    output += '\ttrace_x\ttrace_y\ttrace_z\tsystem'
  else:
    output += '\t' + '\t'.join(['(%i)tx\tty\ttz'%(i+1) for i in range(Nslipsystems[options.lattice])])
output += '\n'

for line in data:
  items = line.split()[:len(labels)]
  if items == []:
    continue
  phi1 = math.radians(float(items[phi1Column]))
  Phi = math.radians(float(items[PhiColumn]))
  phi2 = math.radians(float(items[phi2Column]))
  S = [ sum( [applyEulers(phi1,Phi,phi2,normalize(   slipnormal[options.lattice][slipsystem]))[i]*options.stressnormal[i] for i in range(3)] ) * \
        sum( [applyEulers(phi1,Phi,phi2,normalize(slipdirection[options.lattice][slipsystem]))[i]*options.forcedirection[i] for i in range(3)] ) \
        for slipsystem in range(Nslipsystems[options.lattice]) ]
  output += '\t'.join(items + map(str,S))
  if options.traceplane:
    trace = [crossproduct(options.traceplane,applyEulers(phi1,Phi,phi2,normalize(slipnormal[options.lattice][slipsystem]))) \
               for slipsystem in range(Nslipsystems[options.lattice]) ]
    if options.rank == 0:
      output += '\t' + '\t'.join(map(lambda x:'%f\t%f\t%f'%(x[0],x[1],x[2]),trace))
    elif options.rank > 0:
      SabsSorted = sorted([(abs(S[i]),i) for i in range(len(S))])
      output += '\t' + '\t'.join(map(str,trace[SabsSorted[-options.rank][1]])) + '\t%i'%(1+SabsSorted[-options.rank][1])
      # for t in [normalize(crossproduct(options.traceplane,applyEulers(phi1,Phi,phi2,normalize(slipnormal[options.lattice][i])))) for i in range(12,24)]:
        # print '\t'.join(map(str,t))
        # print '\t'.join(map(lambda x: str(-x),t))
        # print '\t'.join(['0','0','0'])
      # print
  output += '\n'

if filename == []:
  print output
else:
  file = open(filename[0],'w')
  file.write(output)
  file.close()
