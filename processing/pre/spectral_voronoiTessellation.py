#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os,sys,math,string,numpy, damask
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP 


# -----------------------------
class extendedOption(Option):
# -----------------------------
# used for definition of new option parser action 'extend', which enables to take multiple option arguments
# taken from online tutorial http://docs.python.org/library/optparse.html
  
  ACTIONS = Option.ACTIONS + ("extend",)
  STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
  TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
  ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

  def take_action(self, action, dest, opt, value, values, parser):
    if action == "extend":
      lvalue = value.split(",")
      values.ensure_value(dest, []).extend(lvalue)
    else:
      Option.take_action(self, action, dest, opt, value, values, parser)
  



# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
parser = OptionParser(option_class=extendedOption, usage='%prog options [file[s]]', description = """
generates geom file and material_config file using seeds file

""" + string.replace('$Id$','\n','\\n')
)


parser.add_option('-d', '--dimension', dest='dimension', type='float', nargs = 3, \
                                       help='x,y,z dimension of specimen')
parser.add_option('-r', '--resolution', dest='resolution', type='int', nargs = 3, \
                                       help='a,b,c resolution of specimen')
parser.add_option('-o', '--outputName', dest='outputName', type='string', nargs = 1, \
                                       help='Output Name')
                                   
parser.set_defaults(resolution = (0,0,0))
parser.set_defaults(dimension  = (0.0,0.0,0.0)) 
parser.set_defaults(outputName = '')                                    
                                   
(options,filenames) = parser.parse_args()

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    if os.path.splitext(name)[1]=='': name = name+'.seeds'
    if options.outputName=='': options.outputName = name
    if os.path.exists(name):
      files.append({'name':name,\
                    'input':open(name),\
                    'geom':open(os.path.splitext(options.outputName)[0]+'.geom','w+'),\
                    'material.config':open(os.path.splitext(options.outputName)[0]+'_material.config','w+')}) 

# ------------------------------------------ loop over input files ---------------------------------------  
for file in files:
  spatialDim  = 3
  unmapped    = 0

  if file['name'] != 'STDIN': print file['name']
  Favg = numpy.array(([1.0,0.0,0.0],\
                      [0.0,1.0,0.0],\
                      [0.0,0.0,1.0]),'d')

  for lineCount, line in enumerate(file['input']):
    words = line.split()
    if('head' in line): headCount=int(words[0])+1  
    if(lineCount<=headCount):
      if('grains'.lower()==words[0]):
        N_seeds=int(words[1])
        checkGrain=numpy.zeros(N_seeds)
      if 'resolution'==words[0] :
        if options.resolution==(0,0,0):
          resolution=[int(words[2]),int(words[4]),int(words[6])]
        else:
          resolution = [options.resolution[0],options.resolution[1],options.resolution[2]]
       
      if lineCount==headCount:  # all header info there, allocating arrays
        if resolution[2]==1: spatialDim=2
        validDim    = numpy.zeros(spatialDim)
        for i in xrange(spatialDim):
          if(options.dimension[i]>0): validDim[i]= 1
        for i in xrange(spatialDim):  
          if(any(validDim) and options.dimension[i]==0): 
            dimension[i] = max(numpy.array(options.dimension,dtype='float')/resolution) *resolution[i] 
            print 'rescaling invalid dimension '+str(i+1)           
        if(any(validDim)==0):        
          print 'No valid dimension specified, rescaling all dimensions'
          for i in xrange(spatialDim): dimension[i]=1.0/max(resolution)*resolution[i]
        if(all(validDim)): dimension= [options.dimension[0],options.dimension[1],options.dimension[2]] 
        coords = numpy.zeros((spatialDim,N_seeds),'d')
        eulers = numpy.zeros((N_seeds,3),'d')
       
    else:
      Npoints = resolution[0]*resolution[1]*resolution[2]
    
      coords[0:spatialDim,lineCount-headCount] = numpy.array(words[0:spatialDim],'d')
      eulers[lineCount-headCount,0:3] = words[3:6]
      undeformed = numpy.zeros((spatialDim,Npoints),'d').reshape(spatialDim,Npoints) 
      
  for i in xrange(spatialDim):
    for j in xrange(N_seeds):
      if(coords[i][j]>=1.0 or coords[i][j]<0.0):
        print 'WARNING: Seed Coordinate '+ str(coords[i][j]) +' located outside Grid '
  digits = 1+int(math.log10(int(N_seeds)))
  
  file['geom'].write('3 header' +'\n')
  file['geom'].write('resolution  a '+ str(resolution[0])+ '  b '+str(resolution[1])+ '  c '+str(  resolution[2])+'\n')
  file['geom'].write('dimension   x '+ str(dimension[0]) + '  y '+ str(dimension[1])+ '  z '+ str( dimension[2])+'\n')
  file['geom'].write('homogenization  1'+'\n')
  
  file['material.config'].write('<microstructure>'+'\n')

  for i in xrange(N_seeds):
     
        
    file['material.config'].write('[Grain' +str( i+1).zfill(digits)+']'+'\n')
    file['material.config'].write('crystallite 1'+'\n')
    file['material.config'].write('(constituent)  phase 1   texture '+ str(i+1).zfill(digits)+ '   fraction 1.0'+'\n')
  file['material.config'].write('\n'+'<texture>'+'\n')  
  for i in xrange(N_seeds): 
    file['material.config'].write('[Grain'+ str(i+1).zfill(digits)+ ']'+'\n')  
    file['material.config'].write('(gauss)  phi1 '+ str(eulers[i][0])+ '    Phi '+ str(eulers[i][1])+ \
                                        '    Phi2 '+ str(eulers[i][2])+ '   scatter 0.0   fraction 1.0'+'\n')
    
                                                                                 
 
  shift =  dimension[0:spatialDim]/numpy.array( resolution,dtype='float')[0:spatialDim]*0.5  # shift by half of side length to center of element
  for i in xrange(Npoints):
    undeformed[0,i] =  dimension[0]\
                       * float(i                                               % resolution[0])\
                                                                         /float( resolution[0])
    undeformed[1,i] =  dimension[1]\
                       * float(i// resolution[0]                        % resolution[1])\
                                                                         /float( resolution[1])
                                                                      
    if spatialDim==3:
      undeformed[2,i] =  dimension[2]\
                       * float(i// resolution[0]// resolution[1] % resolution[2])\
                                                                         /float( resolution[2])
    undeformed[0:spatialDim,i] += shift
    
    
  indices=damask.core.math.math_nearestNeighborSearch(spatialDim,Favg,numpy.array( dimension,dtype='float'),Npoints,N_seeds,undeformed,coords)
  
  for n in xrange(Npoints):
    file['geom'].write( str(indices[n]//3**spatialDim+1).zfill(digits)+{True:'\n',False:' '}[(n+1)% resolution[0] == 0])
  grainIndex=indices//3**spatialDim+1
  for i in xrange(1,N_seeds+1):
    if i not in grainIndex: unmapped+=1 
  if(unmapped == 0): print 'All Grains Mapped'
  else: print 'Only '+str(N_seeds-unmapped )+' Grains mapped'

  
  file['geom'].close()
  file['material.config'].close()
  file['input'].close()
      