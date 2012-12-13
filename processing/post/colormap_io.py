#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

def write_gsmh(RGB_vector,name):
  colormap = open(str(name) + '.map',"w")
  colormap.write('View.ColorTable = {\n')
  for i in range(len(RGB_vector)-1):
    colormap.write('{'+str((RGB_vector[0][i])*255.0)+','+str((RGB_vector[0][i])*255.0)+','+str((RGB_vector[0][i])*255.0)+'},\n')
  colormap.write('{'+str((RGB_vector[0][-1])*255.0)+','+str((RGB_vector[0][-1])*255.0)+','+str((RGB_vector[0][-1])*255.0)+'}}')
  file.close(colormap)

def write_paraview(RGB_vector,name):
  colormap = open(str(name) + '.xml',"w")
  colormap.write('<ColorMap name = "'+ str(name)+ '" space = "RGB">\n')
  for i in range(len(RGB_vector)):
    colormap.write('<Point x="'+str(i)+'" o="1" r="'+str(RGB_vector[0][i])+'" g="'+str(RGB_vector[1][i])+'" b="'+str(RGB_vector[2][i])+'"/>\n')
  colormap.write('</ColorMap>')
  file.close(colormap)
    
def write_paraview2(RGB_vector,name):
  colormap = open(str(name) + '.xml',"w")
  colormap.write('<ColorMap name = "'+ str(name)+ '" space = "RGB">\n')
  for i in range(len(RGB_vector)/3):
    colormap.write('<Point x="'+str(i)+'" o="1" r="'+str(RGB_vector[i*3])+'" g="'+str(RGB_vector[i*3+1])+'" b="'+str(RGB_vector[i*3+2])+'"/>\n')
  colormap.write('</ColorMap>')
  file.close(colormap)
    
def write_raw(RGB_vector,name):
  colormap = open(str(name) + '.colormap',"w")
  colormap.write('ColorMap name = ' + str(name)+'\n')
  for i in range(len(RGB_vector)):
    colormap.write(str(RGB_vector[0][i])+'\t'+str(RGB_vector[1][i])+'\t'+str(RGB_vector[2][i])+'\n')
  file.close(colormap)

  def read_raw(filename):
    print 'void'
