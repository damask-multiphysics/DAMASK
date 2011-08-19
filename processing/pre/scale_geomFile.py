#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-
#scale the size of a given .geom file (first argument, without extension) in x y z direction by given commandline argument
#(integer in each direction)
#be carefull, last line should be empty (line break after last number)
import sys, time, re
scale = [0 for i in range(3)]
scale[0] = int(sys.argv[2])
scale[1] = int(sys.argv[3])
scale[2] = int(sys.argv[4])

resolution = [3 for i in range(3)]
dimension = [0.0 for i in range(3)]

file_in = open(sys.argv[1]+'.geom','r')
for i in range(3):
  input = file_in.readline()
  if re.match('resolution',input):
    resolution[0] = int(re.findall('\S*',''.join(re.findall('a\s*\S*',input)))[-2])
    resolution[1] = int(re.findall('\S*',''.join(re.findall('b\s*\S*',input)))[-2])
    resolution[2] = int(re.findall('\S*',''.join(re.findall('c\s*\S*',input)))[-2])
  if re.match('dimension',input):
    dimension[0] = float(re.findall('\S*',''.join(re.findall('x\s*\S*',input)))[-2])
    dimension[1] = float(re.findall('\S*',''.join(re.findall('y\s*\S*',input)))[-2])
    dimension[2] = float(re.findall('\S*',''.join(re.findall('z\s*\S*',input)))[-2])
   
output = ''
for x in range(resolution[0]*resolution[1]*resolution[2]):
  a=file_in.readline()
  for i in range(scale[0]):
	  output += a
output_separated = output.splitlines(True)
resolution[0] = resolution[0]*scale[0]    #scale resolution 1 
output = ''
for z in range(resolution[2]):
  for y in range(resolution[1]):
    for i in range(scale[1]):
      output += ''.join(output_separated[z*resolution[0]*resolution[1]+y*resolution[0]:\
                                         z*resolution[0]*resolution[1]+y*resolution[0]+resolution[0]])
resolution[1] = resolution[1]*scale[1]    #scale resolution 2 
output_separated = output.splitlines(True)
output = ''
for z in range(resolution[2]):
  for i in range(scale[2]):
    output += ''.join(output_separated[z*resolution[0]*resolution[1]:\
                                       z*resolution[0]*resolution[1]+resolution[0]*resolution[1]])
resolution[2] = resolution[2]*scale[2]    #scale resolution 2 

file_out = open(sys.argv[1]+'_scaled_'+str(sys.argv[2])+'-'+str(sys.argv[3])+'-'+str(sys.argv[4])+'.geom','w')

file_out.write('resolution a '+str(resolution[0])+' b '+str(resolution[1])+' z '+str(resolution[2])+'\n')
file_out.write('dimension x '+str(dimension[0])+' y '+str(dimension[1])+' z '+str(dimension[2])+'\n')
file_out.write('homogenization 1\n')
file_out.write(output)