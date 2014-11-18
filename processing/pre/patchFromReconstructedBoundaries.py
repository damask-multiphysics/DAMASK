#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import sys,os,math,re,string
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]


try:                                        # check for Python Image Lib
  from PIL import Image,ImageDraw
  ImageCapability = True
except:
  ImageCapability = False

sys.path.append(damask.solver.Marc().libraryPath('../../'))

try:                                        # check for MSC.Mentat Python interface
  from py_mentat import *
  MentatCapability = True
except:
  MentatCapability = False


def outMentat(cmd,locals):
  if cmd[0:3] == '(!)':
    exec(cmd[3:])
  elif cmd[0:3] == '(?)':
    cmd = eval(cmd[3:])
    py_send(cmd)
    if 'log' in locals: locals['log'].append(cmd)
  else:
    py_send(cmd)
    if 'log' in locals: locals['log'].append(cmd)
  return

def outStdout(cmd,locals):
  if cmd[0:3] == '(!)':
    exec(cmd[3:])
  elif cmd[0:3] == '(?)':
    cmd = eval(cmd[3:])
    print cmd
  else:
    print cmd
  return


def output(cmds,locals,dest):
  for cmd in cmds:
    if isinstance(cmd,list):
      output(cmd,locals,dest)
    else:
      {\
      'Mentat': outMentat,\
      'Stdout': outStdout,\
      }[dest](cmd,locals)
  return


def rcbOrientationParser(content,idcolumn):

  grains = []
  myOrientation = [0.0,0.0,0.0]
  for line in content:
    m = re.match(r'\s*(#|$)',line)
    if m: continue                                                # skip comments and blank lines
    for grain in range(2):
      myID = int(line.split()[idcolumn+grain])                    # get grain id
      myOrientation = map(float,line.split())[3*grain:3+3*grain]  # get orientation
      if len(grains) < myID:
        for i in range(myID-len(grains)):                         # extend list to necessary length
          grains.append([0.0,0.0,0.0])
      try:
        grains[myID-1] = myOrientation                            # store Euler angles
      except IndexError:
        message = 'You might not have chosen the correct column for the grain IDs! Please check the "--id" option.'
        print '\033[1;31m'+message+'\033[0m\n'
        raise
      except:
        raise

  return grains

def rcbParser(content,M,size,tolerance,idcolumn,segmentcolumn):   # parser for TSL-OIM reconstructed boundary files

# find bounding box

  boxX = [1.*sys.maxint,-1.*sys.maxint]
  boxY = [1.*sys.maxint,-1.*sys.maxint]
  x = [0.,0.]
  y = [0.,0.]
  for line in content:
    m = re.match(r'\s*(#|$)',line)
    if m: continue                                            # skip comments and blank lines
    try:
      (x[0],y[0],x[1],y[1]) = map(float,line.split())[segmentcolumn:segmentcolumn+4] # get start and end coordinates of each segment.
    except IndexError:
      message = 'You might not have chosen the correct column for the segment end points! Please check the "--segment" option.'
      print '\033[1;31m'+message+'\033[0m\n'
      raise
    except:
      raise
    (x[0],y[0]) = (M[0]*x[0]+M[1]*y[0],M[2]*x[0]+M[3]*y[0])   # apply transformation to coordinates
    (x[1],y[1]) = (M[0]*x[1]+M[1]*y[1],M[2]*x[1]+M[3]*y[1])   # to get rcb --> Euler system
    boxX[0] = min(boxX[0],x[0],x[1])
    boxX[1] = max(boxX[1],x[0],x[1])
    boxY[0] = min(boxY[0],y[0],y[1])
    boxY[1] = max(boxY[1],y[0],y[1])
  dX = boxX[1]-boxX[0]
  dY = boxY[1]-boxY[0]
  
  if size > 0.0: scalePatch = size/dX
  else: scalePatch = 1.0

# read segments
  segment = 0
  connectivityXY = {"0": {"0":[],"%g"%dY:[],},\
                    "%g"%dX: {"0":[],"%g"%dY:[],},}
  connectivityYX = {"0": {"0":[],"%g"%dX:[],},\
                    "%g"%dY: {"0":[],"%g"%dX:[],},}
  grainNeighbors = []
  
  for line in content:
    m = re.match(r'\s*(#|$)',line)
    if m: continue                                # skip comments and blank lines
    (x[0],y[0],x[1],y[1]) = map(float,line.split())[segmentcolumn:segmentcolumn+4] # get start and end coordinates of each segment.
    (x[0],y[0]) = (M[0]*x[0]+M[1]*y[0],M[2]*x[0]+M[3]*y[0])   # apply transformation to coordinates
    (x[1],y[1]) = (M[0]*x[1]+M[1]*y[1],M[2]*x[1]+M[3]*y[1])   # to get rcb --> Euler system

    x[0] -= boxX[0]                                           # make relative to origin of bounding box
    x[1] -= boxX[0]
    y[0] -= boxY[0]
    y[1] -= boxY[0]
    grainNeighbors.append(map(int,line.split()[idcolumn:idcolumn+2])) # remember right and left grain per segment
    for i in range(2):                                      # store segment to both points
      match = False                                         # check whether point is already known (within a small range)
      for posX in connectivityXY.keys():
        if (abs(float(posX)-x[i])<dX*tolerance):
          for posY in connectivityXY[posX].keys():
            if (abs(float(posY)-y[i])<dY*tolerance):
              keyX = posX
              keyY = posY
              match = True
              break
          break
      if (not match):
      # force to boundary if inside tolerance to it
        if (abs(x[i])<dX*tolerance):
          x[i] = 0
        if (abs(dX-x[i])<dX*tolerance):
          x[i] = dX
        if (abs(y[i])<dY*tolerance):
          y[i] = 0
        if (abs(dY-y[i])<dY*tolerance):
          y[i] = dY
        keyX = "%g"%x[i]
        keyY = "%g"%y[i]
        if keyX not in connectivityXY:                      # create new hash entry for so far unknown point
          connectivityXY[keyX] = {}
        if keyY not in connectivityXY[keyX]:                # create new hash entry for so far unknown point
          connectivityXY[keyX][keyY] = []
        if keyY not in connectivityYX:                      # create new hash entry for so far unknown point
          connectivityYX[keyY] = {}
        if keyX not in connectivityYX[keyY]:                # create new hash entry for so far unknown point
          connectivityYX[keyY][keyX] = []
      connectivityXY[keyX][keyY].append(segment)
      connectivityYX[keyY][keyX].append(segment)
    segment += 1
      

# top border
  keyId = "0"
  boundary = connectivityYX[keyId].keys()
  boundary.sort(key=float)
  for indexBdy in range(len(boundary)-1):
    connectivityXY[boundary[indexBdy]][keyId].append(segment)
    connectivityXY[boundary[indexBdy+1]][keyId].append(segment)
    connectivityYX[keyId][boundary[indexBdy]].append(segment)
    connectivityYX[keyId][boundary[indexBdy+1]].append(segment)
    segment += 1

# right border
  keyId = "%g"%(boxX[1]-boxX[0])
  boundary = connectivityXY[keyId].keys()
  boundary.sort(key=float)
  for indexBdy in range(len(boundary)-1):
    connectivityYX[boundary[indexBdy]][keyId].append(segment)
    connectivityYX[boundary[indexBdy+1]][keyId].append(segment)
    connectivityXY[keyId][boundary[indexBdy]].append(segment)
    connectivityXY[keyId][boundary[indexBdy+1]].append(segment)
    segment += 1

# bottom border
  keyId = "%g"%(boxY[1]-boxY[0])
  boundary = connectivityYX[keyId].keys()
  boundary.sort(key=float,reverse=True)
  for indexBdy in range(len(boundary)-1):
    connectivityXY[boundary[indexBdy]][keyId].append(segment)
    connectivityXY[boundary[indexBdy+1]][keyId].append(segment)
    connectivityYX[keyId][boundary[indexBdy]].append(segment)
    connectivityYX[keyId][boundary[indexBdy+1]].append(segment)
    segment += 1

# left border
  keyId = "0"
  boundary = connectivityXY[keyId].keys()
  boundary.sort(key=float,reverse=True)
  for indexBdy in range(len(boundary)-1):
    connectivityYX[boundary[indexBdy]][keyId].append(segment)
    connectivityYX[boundary[indexBdy+1]][keyId].append(segment)
    connectivityXY[keyId][boundary[indexBdy]].append(segment)
    connectivityXY[keyId][boundary[indexBdy+1]].append(segment)
    segment += 1


  allkeysX = connectivityXY.keys()
  allkeysX.sort()
  points = []
  segments = [[] for i in range(segment)]
  pointId = 0
  for keyX in allkeysX:
    allkeysY = connectivityXY[keyX].keys()
    allkeysY.sort()
    for keyY in allkeysY:
      points.append({'coords': [float(keyX)*scalePatch,float(keyY)*scalePatch], 'segments': connectivityXY[keyX][keyY]})
      for segment in connectivityXY[keyX][keyY]:
        if (segments[segment] == None):
          segments[segment] = pointId
        else:
          segments[segment].append(pointId)
      pointId += 1

  crappyData = False
  for pointId,point in enumerate(points):
    if len(point['segments']) < 2:                            # point marks a dead end!
      print "Dead end at segment %i (%f,%f)"%(1+point['segments'][0],boxX[0]+point['coords'][0]/scalePatch,boxY[0]+point['coords'][1]/scalePatch,)
      crappyData = True

  grains = {'draw': [], 'legs': []}

  if not crappyData:

    for pointId,point in enumerate(points):
      while point['segments']:
        myStart = pointId
        grainDraw = [points[myStart]['coords']]
        innerAngleSum = 0.0
        myWalk = point['segments'].pop()
        grainLegs = [myWalk]
        if segments[myWalk][0] == myStart:
          myEnd = segments[myWalk][1]
        else:
          myEnd = segments[myWalk][0]

        while (myEnd != pointId):
          myV = [points[myEnd]['coords'][0]-points[myStart]['coords'][0],\
            points[myEnd]['coords'][1]-points[myStart]['coords'][1]]
          myLen = math.sqrt(myV[0]**2+myV[1]**2)
          best = {'product': -2.0, 'peek': -1, 'len': -1, 'point': -1}
          for peek in points[myEnd]['segments']:                                                    # trying in turn all segments emanating from current endPoint
            if peek == myWalk:
              continue
            peekEnd = segments[peek][1] if segments[peek][0] == myEnd else segments[peek][0]
            peekV = [points[myEnd]['coords'][0]-points[peekEnd]['coords'][0],
                     points[myEnd]['coords'][1]-points[peekEnd]['coords'][1]]
            peekLen = math.sqrt(peekV[0]**2+peekV[1]**2)
            crossproduct = (myV[0]*peekV[1]-myV[1]*peekV[0])/myLen/peekLen
            dotproduct = (myV[0]*peekV[0]+myV[1]*peekV[1])/myLen/peekLen
            if crossproduct*(dotproduct+1.0) >= best['product']:
              best['product'] = crossproduct*(dotproduct+1.0)
              best['peek'] = peek
              best['point'] = peekEnd
          innerAngleSum += best['product']
          myWalk = best['peek']
          myStart = myEnd
          myEnd = best['point']
          if myWalk in points[myStart]['segments']:
            points[myStart]['segments'].remove(myWalk)
          else:
            sys.stderr.write(str(myWalk)+' not in segments of point '+str(myStart)+'\n')
          grainDraw.append(points[myStart]['coords'])
          grainLegs.append(myWalk)
        if innerAngleSum > 0.0:
          grains['draw'].append(grainDraw)
          grains['legs'].append(grainLegs)
        else:
          grains['box'] = grainLegs

  
# build overall data structure

  rcData = {'dimension':[dX,dY], 'point': [],'segment': [], 'grain': [], 'grainMapping': []}
  print "  dimension %g x %g"%(dX,dY)

  for point in points:
    rcData['point'].append(point['coords'])
  print "  found %i points"%(len(rcData['point']))

  for segment in segments:  
    rcData['segment'].append(segment)
  print "  built %i segments"%(len(rcData['segment']))

  for legs in grains['legs']:                                                             # loop over grains
    rcData['grain'].append(legs)                                                          # store list of boundary segments
    myNeighbors = {}
    for leg in legs:                                                                      # test each boundary segment
      if leg < len(grainNeighbors):                                                       # a valid segment index?
        for side in range(2):                                                             # look at both sides of the segment
          if grainNeighbors[leg][side] in myNeighbors:                                    # count occurrence of grain IDs
            myNeighbors[grainNeighbors[leg][side]] += 1
          else:
            myNeighbors[grainNeighbors[leg][side]] = 1
    if myNeighbors:                                                                       # do I have any neighbors (i.e., non-bounding box segment)
      candidateGrains = sorted(myNeighbors.iteritems(), key=lambda (k,v): (v,k), reverse=True)  # sort grain counting
      if candidateGrains[0][0] not in rcData['grainMapping']:                                   # most frequent one not yet seen?
        rcData['grainMapping'].append(candidateGrains[0][0])                                    # must be me then
      else:
        rcData['grainMapping'].append(candidateGrains[1][0])                                    # special case of bi-crystal situation...
      
  print "  found %i grains\n"%(len(rcData['grain']))

  rcData['box'] = grains['box'] if 'box' in grains else []
  
  return rcData


def init():
    return ["*new_model yes",
      "*select_clear",
      "*reset",
      "*set_nodes off",
      "*elements_solid",
      "*show_view 4",
      "*reset_view",
      "*view_perspective",
      "*redraw",
      ]


def sample(size,aspect,n,xmargin,ymargin):

  cmds = [\
# gauge
    "*add_points %f %f %f"%(-size*(0.5+xmargin), size*(0.5*aspect+ymargin),0),
    "*add_points %f %f %f"%( size*(0.5+xmargin), size*(0.5*aspect+ymargin),0),
    "*add_points %f %f %f"%( size*(0.5+xmargin),-size*(0.5*aspect+ymargin),0),
    "*add_points %f %f %f"%(-size*(0.5+xmargin),-size*(0.5*aspect+ymargin),0),
    "*set_curve_type line",
    "*add_curves %i %i"%(1,2),
    "*add_curves %i %i"%(3,4),
    "*set_curve_div_type_fix_ndiv",
    "*set_curve_div_num %i"%n,
    "*apply_curve_divisions",
    "1 2 #",
    "*add_curves %i %i"%(2,3),  # right side
    "*add_curves %i %i"%(4,1),  # left side
    "*set_curve_div_type_fix_ndiv",
    "*set_curve_div_num %i"%n,
    "*apply_curve_divisions",
    "3 4 #",
    ]

  return cmds
  

def patch(a,n,mesh,rcData):
  cmds = []
  for l in range(len(rcData['point'])):           # generate all points
    cmds.append("*add_points %f %f %f"%(rcData['point'][l][0]-a/2.0,rcData['point'][l][1]-a/rcData['dimension'][0]*rcData['dimension'][1]/2.0,0))

  cmds.append(["*set_curve_type line",
         "*set_curve_div_type_fix_ndiv",
        ])
  for m in range(len(rcData['segment'])):           # generate all curves and subdivide them for overall balanced piece length
    start = rcData['segment'][m][0]
    end   = rcData['segment'][m][1]
    cmds.append([\
      "*add_curves %i %i" %(start+rcData['offsetPoints'],
                  end  +rcData['offsetPoints']),
      "*set_curve_div_num %i"%(max(1,round(math.sqrt((rcData['point'][start][0]-rcData['point'][end][0])**2+\
                               (rcData['point'][start][1]-rcData['point'][end][1])**2)/a*n))),
      "*apply_curve_divisions",
      "%i #"%(m+rcData['offsetSegments']),
      ])

  grain = 0
  cmds.append('(!)locals["last"] = py_get_int("nelements()")')
  for g in rcData['grain']:
    cmds.append([\
      '(!)locals["first"] = locals["last"]+1',
      "*%s "%mesh+" ".join([str(rcData['offsetSegments']+x) for x in g])+" #",
      '(!)locals["last"] = py_get_int("nelements()")',
      "*select_elements",
      '(?)"%i to %i #"%(locals["first"],locals["last"])',
      "*store_elements grain_%i"%rcData['grainMapping'][grain],
      "all_selected",
      "*select_clear",
      ])
    grain += 1

  return cmds


def gage(mesh,rcData):

  return([\
    "*%s "%mesh + 
    " ".join([str(x) for x in range(1,rcData['offsetSegments'])]) +
    " " +
    " ".join([str(rcData['offsetSegments']+x)for x in rcData['box']]) +
    " #",
    "*select_reset",
    "*select_clear",
    "*select_elements",
    "all_existing",
    "*select_mode_except",
    ['grain_%i'%rcData['grainMapping'][i] for i in range(len(rcData['grain']))],
    "#",
    "*store_elements matrix",
    "all_selected",
    "*select_mode_invert",
    "*select_elements",
    "all_existing",
    "*store_elements _grains",
    "all_selected",
    "*select_clear",
    "*select_reset",
    ])


def expand3D(thickness,steps):
  return([\
    "*set_expand_translation z %f"%(thickness/steps),
    "*set_expand_repetitions %i"%steps,
    "*expand_elements",
    "all_existing",
    ])


def initial_conditions(grainNumber,grainMapping):
  cmds = [\
    "*new_icond",
    "*icond_name _temperature",
    "*icond_type state_variable",
    "*icond_param_value state_var_id 1",
    "*icond_dof_value var 300",
    "*add_icond_elements",
    "all_existing",
    "*new_icond",
    "*icond_name _homogenization",
    "*icond_type state_variable",
    "*icond_param_value state_var_id 2",
    "*icond_dof_value var 1",
    "*add_icond_elements",
    "all_existing",
    ]
  
  for grain in range(grainNumber):
    cmds.append([\
    "*new_icond",
    "*icond_name grain_%i"%grainMapping[grain],
    "*icond_type state_variable",
    "*icond_param_value state_var_id 3",
    "*icond_dof_value var %i"%(grain+1),
    "*add_icond_elements",
    "grain_%i"%grainMapping[grain],
    "",
          ])
  cmds.append([\
    "*new_icond",
    "*icond_name rim",
    "*icond_type state_variable",
    "*icond_param_value state_var_id 3",
    "*icond_dof_value var %i"%(grainNumber+1),
    "*add_icond_elements",
    "matrix",
    ])
  return cmds


def boundary_conditions(rate,thickness, size,aspect,xmargin,ymargin):

  inner = (1 - 1.0e-4) * size*(0.5+xmargin)
  outer = (1 + 1.0e-4) * size*(0.5+xmargin)
  lower = (1 - 1.0e-4) * size*(0.5*aspect+ymargin)
  upper = (1 + 1.0e-4) * size*(0.5*aspect+ymargin)
  
  return [\
  "*new_md_table 1 1",
  "*table_name linear",
  "*set_md_table_type 1 time",
  "*table_add",
  "0 0",
  "1 1",
  "*select_method_box",
  "*new_apply",
  "*apply_name pull_bottom",
  "*apply_type fixed_displacement",
  "*apply_dof y",
  "*apply_dof_value y %f"%(-rate*(lower+upper)/2.0),
  "*apply_dof_table y linear",
  "*select_clear_nodes",
  "*select_nodes",
  "%f %f"%(-outer,outer),
  "%f %f"%(-upper,-lower),
  "%f %f"%(-.0001*thickness,1.0001*thickness),
  "*add_apply_nodes",
  "all_selected",
  "*new_apply",
  "*apply_name pull_top",
  "*apply_type fixed_displacement",
  "*apply_dof y",
  "*apply_dof_value y %f"%(rate*(lower+upper)/2.0),
  "*apply_dof_table y linear",
  "*select_clear_nodes",
  "*select_nodes",
  "%f %f"%(-outer,outer),
  "%f %f"%(lower,upper),
  "%f %f"%(-.0001*thickness,1.0001*thickness),
  "*add_apply_nodes",
  "all_selected",
  "*new_apply",
  "*apply_name fix_x",
  "*apply_type fixed_displacement",
  "*apply_dof x",
  "*apply_dof_value x 0",
  "*select_clear_nodes",
  "*select_nodes",
  "%f %f"%(-outer,-inner),
  "%f %f"%(lower,upper),
  "%f %f"%(-.0001*thickness,.0001*thickness),
  "%f %f"%(-outer,-inner),
  "%f %f"%(lower,upper),
  "%f %f"%(0.9999*thickness,1.0001*thickness),
  "%f %f"%(-outer,-inner),
  "%f %f"%(-upper,-lower),
  "%f %f"%(-.0001*thickness,.0001*thickness),
  "%f %f"%(-outer,-inner),
  "%f %f"%(-upper,-lower),
  "%f %f"%(0.9999*thickness,1.0001*thickness),
  "*add_apply_nodes",
  "all_selected",
  "*new_apply",
  "*apply_name fix_z",
  "*apply_type fixed_displacement",
  "*apply_dof z",
  "*apply_dof_value z 0",
  "*select_clear_nodes",
  "*select_nodes",
  "%f %f"%(-outer,-inner),
  "%f %f"%(lower,upper),
  "%f %f"%(-.0001*thickness,.0001*thickness),
  "%f %f"%(-outer,-inner),
  "%f %f"%(-upper,-lower),
  "%f %f"%(-.0001*thickness,.0001*thickness),
  "%f %f"%(inner,outer),
  "%f %f"%(lower,upper),
  "%f %f"%(-.0001*thickness,.0001*thickness),
  "%f %f"%(inner,outer),
  "%f %f"%(-upper,-lower),
  "%f %f"%(-.0001*thickness,.0001*thickness),
  "*add_apply_nodes",
  "all_selected",
  "*select_clear",
  "*select_reset",
  ]

def materials():
  return [\
  "*new_material",
  "*material_name patch",
  "*material_type mechanical:hypoelastic",
  "*material_option hypoelastic:method:hypela2",
  "*material_option hypoelastic:pass:def_rot",
  "*add_material_elements",
  "all_existing",
  ]
  

def loadcase(time,incs,Ftol): 
  return [\
  "*new_loadcase",
  "*loadcase_name puller",
  "*loadcase_type static",
  "*loadcase_value time",
  "%g"%time,
  "*loadcase_value nsteps",
  "%i"%incs,
  "*loadcase_value maxrec",
  "20",
  "*loadcase_value ntime_cuts",
  "30",
  "*loadcase_value force",
  "%g"%Ftol,
  ]


def job(grainNumber,grainMapping,twoD):  
  return [\
  "*new_job",
  "*job_name pull",
  "*job_class mechanical",
  "*add_job_loadcases puller",
  "*add_job_iconds homogenization",
  ["*add_job_iconds grain_%i"%i for i in grainMapping[:grainNumber]],
  "*add_job_iconds rim",
  "*job_option dimen:%s                  | analysis dimension"%({True:'two  ',False:'three'}[twoD]),
  "*job_option strain:large                 | finite strains",
  "*job_option large_strn_proc:upd_lagrange | updated Lagrange framework",
  "*job_option plas_proc:multiplicative     | multiplicative decomp of F",
  "*job_option solver_nonsym:on             | nonsymmetrical solution",
  "*job_option solver:mfront_sparse         | multi-frontal sparse",
  "*job_param stef_boltz 5.670400e-8",
  "*job_param univ_gas_const 8.314472",
  "*job_param planck_radiation_2 1.4387752e-2",
  "*job_param speed_light_vacuum 299792458",
#  "*job_usersub_file /san/%s/FEM/DAMASK/code/mpie_cpfem_marc2010.f90  | subroutine definition"%(pwd.getpwuid(os.geteuid())[0].rpartition("\\")[2]),
  "*job_option user_source:compile_save",
  ]

# "*job_option large:on                  | large displacement",
# "*job_option plasticity:l_strn_mn_add  | large strain additive",
# "*job_option cdilatation:on            | constant dilatation",
# "*job_option update:on                 | updated lagrange procedure",
# "*job_option finite:on                 | large strains",
# "*job_option restart_mode:write        | enable restarting",


def postprocess():
  return [\
  "*add_post_tensor stress",
  "*add_post_tensor strain",
  "*add_post_var von_mises",
  "",
  ]



def cleanUp(a):
  return [\
  "*remove_curves",
  "all_existing",
  "*remove_points",
  "all_existing",
  "*set_sweep_tolerance %f"%(1e-5*a),
  "*sweep_all",
  "*renumber_all",
  ]


# -------------------------
def image(name,imgsize,marginX,marginY,rcData):
# -------------------------

  dX = max([coords[0] for coords in rcData['point']])
  dY = max([coords[1] for coords in rcData['point']])
  offsetX = imgsize*marginX
  offsetY = imgsize*marginY
  sizeX = int(imgsize*(1    +2*marginX))
  sizeY = int(imgsize*(dY/dX+2*marginY))

  scaleImg = imgsize/dX            # rescale from max x coord

  img = Image.new("RGB",(sizeX,sizeY),(255,255,255))
  draw = ImageDraw.Draw(img)

  for id,point in enumerate(rcData['point']):
    draw.text([offsetX+point[0]*scaleImg,sizeY-(offsetY+point[1]*scaleImg)],"%i"%id,fill=(0,0,0))

  for id,vertex in enumerate(rcData['segment']):
      start = rcData['point'][vertex[0]]
      end   = rcData['point'][vertex[1]]
      draw.text([offsetX+(start[0]+end[0])/2.0*scaleImg,sizeY-(offsetY+(start[1]+end[1])/2.0*scaleImg)],"%i"%id,fill=(0,0,128))
      draw.line([offsetX+start[0]*scaleImg,sizeY-(offsetY+start[1]*scaleImg),
                 offsetX+  end[0]*scaleImg,sizeY-(offsetY+  end[1]*scaleImg)],width=1,fill=(128,128,128))

  for id,segment in enumerate(rcData['box']):
        start = rcData['point'][rcData['segment'][segment][0]]
        end   = rcData['point'][rcData['segment'][segment][1]]
        draw.line([offsetX+start[0]*scaleImg,sizeY-(offsetY+start[1]*scaleImg),
                   offsetX+  end[0]*scaleImg,sizeY-(offsetY+  end[1]*scaleImg)],width=3,fill=(128,128*(id%2),0))

  for grain,origGrain in enumerate(rcData['grainMapping']):
    center = [0.0,0.0]
    for segment in rcData['grain'][grain]:                    # loop thru segments around grain
      for point in rcData['segment'][segment]:                # take start and end points
        center[0] += rcData['point'][point][0]                # build vector sum
        center[1] += rcData['point'][point][1]

    center[0] /= len(rcData['grain'][grain])*2                 # normalize by two times segment count, i.e. point count
    center[1] /= len(rcData['grain'][grain])*2

    draw.text([offsetX+center[0]*scaleImg,sizeY-(offsetY+center[1]*scaleImg)],'%i -> %i'%(grain,origGrain),fill=(128,32,32))

  img.save(name+'.png',"PNG")

# -------------------------
def inside(x,y,points):                                                     # tests whether point(x,y) is within polygon described by points
# -------------------------
 inside = False
 npoints=len(points)
 (x1,y1) = points[npoints-1]                                                # start with last point of points
 startover = (y1 >= y)                                                      # am I above testpoint?
 for i in range(npoints):                                                   # loop through all points
   (x2,y2) = points[i]                                                      # next point
   endover = (y2 >= y)                                                      # am I above testpoint?
   if (startover != endover):                                               # one above one below testpoint?
     if((y2 - y)*(x2 - x1) <= (y2 - y1)*(x2 - x)):                          # check for intersection
       if (endover):
         inside = not inside                                                # found intersection
     else:
       if (not endover):
         inside = not inside                                                # found intersection
   startover = endover                                                      # make second point first point
   (x1,y1) = (x2,y2)
   
 return inside
  
# -------------------------
def fftbuild(rcData,height,xframe,yframe,resolution,extrusion):            # build array of grain numbers
# -------------------------
 maxX = -1.*sys.maxint
 maxY = -1.*sys.maxint
 for line in rcData['point']:                                               # find data range
   (x,y) = line
   maxX = max(maxX, x)
   maxY = max(maxY, y)
 xsize = maxX+2*xframe                                                      # add framsize
 ysize = maxY+2*yframe
 xres = int(round(resolution/2.0)*2)                                        # use only even resolution
 yres = int(round(xres/xsize*ysize/2.0)*2)                                  # calculate other resolutions
 zres = extrusion
 zsize = extrusion*min([xsize/xres,ysize/yres])
 
 fftdata = {'fftpoints':[], \
            'resolution':(xres,yres,zres), \
            'dimension':(xsize,ysize,zsize)}
 
 frameindex=len(rcData['grain'])+1                                          # calculate frame index as largest grain index plus one 
 dx = xsize/(xres+1)                                                        # calculate step sizes
 dy = ysize/(yres+1)

 grainpoints = []
 for segments in rcData['grain']:                                           # get segments of each grain
   points = {}
   for i,segment in enumerate(segments[:-1]):                               # loop thru segments except last (s=[start,end])
     points[rcData['segment'][segment][0]] = i                              # assign segment index to start point
     points[rcData['segment'][segment][1]] = i                              # assigne segment index to endpoint
   for i in range(2):                                                       # check points of last segment
     if points[rcData['segment'][segments[-1]][i]] != 0:                    # not on first segment
       points[rcData['segment'][segments[-1]][i]] = len(segments)-1         # assign segment index to last point
   
   grainpoints.append([])                                                   # start out blank for current grain
   for p in sorted(points, key=points.get):                                 # loop thru set of sorted points
     grainpoints[-1].append([rcData['point'][p][0],rcData['point'][p][1]])  # append x,y of point

 bestGuess = 0                                                              # assume grain 0 as best guess
 for i in range(int(xres*yres)):                                            # walk through all points in xy plane
   xtest = -xframe+((i%xres)+0.5)*dx                                        # calculate coordinates
   ytest = -yframe+(int(i/xres)+0.5)*dy
   if(xtest < 0 or xtest > maxX):                                           # check wether part of frame
     if( ytest < 0 or ytest > maxY):                                        # part of edges
       fftdata['fftpoints'].append(frameindex+2)                            # append frameindex to result array
     else:                                                                  # part of xframe
       fftdata['fftpoints'].append(frameindex)                              # append frameindex to result array
   elif( ytest < 0 or ytest > maxY):                                        # part of yframe
       fftdata['fftpoints'].append(frameindex+1)                            # append frameindex to result array
   else:
      if inside(xtest,ytest,grainpoints[bestGuess]):                        # check best guess first
        fftdata['fftpoints'].append(bestGuess+1)
      else:                                                                 # no success
        for g in range(len(grainpoints)):                                   # test all
          if inside(xtest,ytest,grainpoints[g]):
            fftdata['fftpoints'].append(g+1)
            bestGuess = g
            break
        
 return fftdata  


# ----------------------- MAIN -------------------------------
parser = OptionParser(option_class=damask.extendableOption, usage='%prog [options] datafile[s]', description = """
Produce image, spectral geometry description, and (auto) Mentat procedure from TSL/OIM
reconstructed boundary file

""", version = scriptID)

parser.add_option("-o", "--output", action='extend', dest='output', \
        help="types of output [image, mentat, procedure, spectral]")
parser.add_option("-p", "--port", type="int",\
        dest="port",\
        help="Mentat connection port [%default]")
parser.add_option("-2", "--twodimensional", action="store_true",\
        dest="twoD",\
        help="twodimensional model [%default]")
parser.add_option("-s","--patchsize", type="float",\
        dest="size",\
        help="height of patch [%default]")
parser.add_option("-e", "--strain", type="float",\
        dest="strain",\
        help="final strain to reach in simulation [%default]")
parser.add_option("--rate", type="float",\
        dest="strainrate",\
        help="(engineering) strain rate to simulate [%default]")
parser.add_option("-N", "--increments", type="int",\
        dest="increments",\
        help="number of increments to take [%default]")
parser.add_option("-t", "--tolerance", type="float",\
        dest="tolerance",\
        help="relative tolerance of pixel positions to be swept [%default]")
parser.add_option("-m", "--mesh", choices=['dt_planar_trimesh','af_planar_trimesh','af_planar_quadmesh'],\
        dest="mesh",\
        help="algorithm and element type for automeshing (dt_planar_trimesh, af_planar_trimesh, af_planar_quadmesh) [%default]")
parser.add_option("-x", "--xmargin", type="float",\
        dest="xmargin",\
        help="margin in x in units of patch size [%default]")
parser.add_option("-y", "--ymargin", type="float",\
        dest="ymargin",\
        help="margin in y in units of patch size [%default]")
parser.add_option("-r", "--resolution", type="int",\
        dest="resolution",\
        help="number of Fourier points/Finite Elements across patch size + x_margin [%default]")
parser.add_option("-z", "--extrusion", type="int",\
        dest="extrusion",\
        help="number of repetitions in z-direction [%default]")
parser.add_option("-i", "--imagesize", type="int",\
        dest="imgsize",\
        help="size of PNG image [%default]")
parser.add_option("-M", "--coordtransformation", type="float", nargs=4, \
        dest="M",\
        help="2x2 transformation from rcb to Euler coords [%default]")
parser.add_option("--scatter", type="float",\
        dest="scatter",\
        help="orientation scatter %default")
parser.add_option("--segment", type="int",\
        dest="segmentcolumn",\
        help="column holding the first entry for the segment end points in the rcb file [%default]")
parser.add_option("--id", type="int",\
        dest="idcolumn",\
        help="column holding the right hand grain ID in the rcb file [%default]")

parser.set_defaults(output = [])
parser.set_defaults(size = 1.0)
parser.set_defaults(port = 40007)
parser.set_defaults(xmargin = 0.0)
parser.set_defaults(ymargin = 0.0)
parser.set_defaults(resolution = 64)
parser.set_defaults(extrusion = 2)
parser.set_defaults(imgsize = 512)
parser.set_defaults(M = [0.0,1.0,1.0,0.0])  # M_11, M_12, M_21, M_22.  x,y in RCB is y,x of Eulers!!
parser.set_defaults(tolerance = 1.0e-3)
parser.set_defaults(scatter = 0.0)
parser.set_defaults(strain = 0.2)
parser.set_defaults(strainrate = 1.0e-3)
parser.set_defaults(increments = 200)
parser.set_defaults(mesh = 'dt_planar_trimesh')
parser.set_defaults(twoD = False)
parser.set_defaults(segmentcolumn = 9)
parser.set_defaults(idcolumn = 13)

(options, args) = parser.parse_args()

print '\033[1m'+scriptName+'\033[0m\n'
if not len(args):
  parser.error('no boundary file specified')

try:
  boundaryFile = open(args[0])
  boundarySegments = boundaryFile.readlines()
  boundaryFile.close()
except:
  print 'unable to read boundary file "%s"'%args[0]
  raise

options.output = [s.lower() for s in options.output]                        # lower case
options.idcolumn -= 1                                                       # python indexing starts with 0
options.segmentcolumn -= 1                                                  # python indexing starts with 0

myName = os.path.splitext(args[0])[0]
print "%s\n"%myName

orientationData = rcbOrientationParser(boundarySegments,options.idcolumn)
rcData = rcbParser(boundarySegments,options.M,options.size,options.tolerance,options.idcolumn,options.segmentcolumn)

# ----- write image -----

if 'image' in options.output and options.imgsize > 0:
  if ImageCapability:
    image(myName,options.imgsize,options.xmargin,options.ymargin,rcData)
  else:
    print '...no image drawing possible (PIL missing)...'

# ----- write spectral geom -----

if 'spectral' in options.output:
  fftdata = fftbuild(rcData, options.size, options.xmargin, options.ymargin, options.resolution, options.extrusion)
  
  geomFile = open(myName+'_'+str(int(fftdata['resolution'][0]))+'.geom','w')  # open geom file for writing
  geomFile.write('3\theader\n')                                               # write header info
  geomFile.write('resolution a %i b %i c %i\n'%(fftdata['resolution']))       # resolution
  geomFile.write('dimension x %f y %f z %f\n'%(fftdata['dimension']))         # size
  geomFile.write('homogenization 1\n')                                        # homogenization
  for z in xrange(fftdata['resolution'][2]):                                  # z repetions
    for y in xrange(fftdata['resolution'][1]):                                # each x-row separately
      geomFile.write('\t'.join(map(str,fftdata['fftpoints'][ y   *fftdata['resolution'][0]:
                                                            (y+1)*fftdata['resolution'][0]]))+'\n')  # grain indexes, x-row per line
  geomFile.close()                                                            # close geom file
  
  print('assigned %i out of %i (2D) Fourier points.'%(len(fftdata['fftpoints']), int(fftdata['resolution'][0])*int(fftdata['resolution'][1])))
  

# ----- write Mentat procedure -----

if 'mentat' in options.output:
  if MentatCapability:
    rcData['offsetPoints']   = 1+4                # gage definition generates 4 points
    rcData['offsetSegments'] = 1+4                # gage definition generates 4 segments
    
    cmds = [\
      init(),
      sample(options.size,rcData['dimension'][1]/rcData['dimension'][0],12,options.xmargin,options.ymargin),
      patch(options.size,options.resolution,options.mesh,rcData),
      gage(options.mesh,rcData),
      ]
    
    if not options.twoD:
      cmds += [expand3D(options.size*(1.0+2.0*options.xmargin)/options.resolution*options.extrusion,options.extrusion),]
    
    cmds += [\
      cleanUp(options.size),
      materials(),
      initial_conditions(len(rcData['grain']),rcData['grainMapping']),
      boundary_conditions(options.strainrate,options.size*(1.0+2.0*options.xmargin)/options.resolution*options.extrusion,\
                          options.size,rcData['dimension'][1]/rcData['dimension'][0],options.xmargin,options.ymargin),
      loadcase(options.strain/options.strainrate,options.increments,0.01),
      job(len(rcData['grain']),rcData['grainMapping'],options.twoD),
      postprocess(),
      ["*identify_sets","*regen","*fill_view","*save_as_model %s yes"%(myName)],
    ]
    
    outputLocals = {'log':[]}
    if (options.port != None):
      py_connect('',options.port)
      try: 
        output(cmds,outputLocals,'Mentat')
      finally:
        py_disconnect()
      if 'procedure' in options.output:
        output(outputLocals['log'],outputLocals,'Stdout')
  else:
    print '...no interaction with Mentat possible...'
    

# ----- write config data to file -----

if 'mentat' in options.output or 'spectral' in options.output:
  output = ''
  output += '\n\n<homogenization>\n' + \
            '\n[SX]\n' + \
            'type\tisostrain\n' + \
            'Ngrains\t1\n' + \
            '\n\n<microstructure>\n'
  
  for i,grain in enumerate(rcData['grainMapping']):
    output += '\n[grain %i]\n'%grain + \
              'crystallite\t1\n' + \
              '(constituent)\tphase 1\ttexture %i\tfraction 1.0\n'%(i+1)
  if (options.xmargin > 0.0):
      output += '\n[x-margin]\n' + \
                'crystallite\t1\n' + \
                '(constituent)\tphase 2\ttexture %i\tfraction 1.0\n'%(len(rcData['grainMapping'])+1)
  if (options.ymargin > 0.0):
      output += '\n[y-margin]\n' + \
                'crystallite\t1\n' + \
                '(constituent)\tphase 2\ttexture %i\tfraction 1.0\n'%(len(rcData['grainMapping'])+1)
  if (options.xmargin > 0.0 and options.ymargin > 0.0):
      output += '\n[margin edge]\n' + \
                'crystallite\t1\n' + \
                '(constituent)\tphase 2\ttexture %i\tfraction 1.0\n'%(len(rcData['grainMapping'])+1)
  
  output += '\n\n<crystallite>\n' + \
            '\n[fillMeIn]\n' + \
            '\n\n<phase>\n' + \
            '\n[patch]\n'
  if (options.xmargin > 0.0 or options.ymargin > 0.0):
    output += '\n[margin]\n'
  
  output += '\n\n<texture>\n\n'
  for grain in rcData['grainMapping']:
    output += '\n[grain %i]\n'%grain + \
              '(gauss)\tphi1\t%f\tphi\t%f\tphi2\t%f\tscatter\t%f\tfraction\t1.0\n'\
              %(math.degrees(orientationData[grain-1][0]),math.degrees(orientationData[grain-1][1]),math.degrees(orientationData[grain-1][2]),options.scatter)
  if (options.xmargin > 0.0 or options.ymargin > 0.0):
    output += '\n[margin]\n' + \
              '(random)\t\tscatter\t0.0\tfraction\t1.0\n'
  
  configFile = open(myName+'.config','w')
  configFile.write(output)
  configFile.close()
