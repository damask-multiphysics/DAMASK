#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

# This script is used to generate colormaps for gmsh (http://geuz.org/gmsh/)
# The script writes 360 files. Each file contains one colormap.
# More information on the used colors space can be found at http://en.wikipedia.org/wiki/HSL_and_HSV
# written by M. Diehl, m.diehl@mpie.de

import math

print '******************************************************************************'
print '                   Write colormaps for gmsh'
print ''
print 'Suitable for datasets that have only positive/negative values'
print 'The colors are described using the HSL model.'
print 'Each of the 360 generated colormaps uses one value of (H)ue.'
print 'The colormaps runs at constant H from given (L)ightness and (S)aturation'
print 'to given L and S'
print 'L is distribute linearly, S changes as the square root of a linear list.'
print 'Suitable values: L_start = S_start = 1 , L_end = S_end = 0.'
print '******************************************************************************'
print ''
startL = float(raw_input('Please enter start value for (L)ightness: '))
endL = float(raw_input('Please enter end value for L: '))
startS = float(raw_input('Please enter start value for (S)aturation: '))
endS = float(raw_input('Please enter end value for S: '))
steps = int(raw_input('Please enter steps/resolution: '))
for h in range(0,360):
	colormap = open('colormap_' + str(h).zfill(3) + '.map',"w")
	colormap.write('View.ColorTable = {\n')
	for i in range(0,steps):
		h_strich = h/60.0
		if(h_strich>6.0):
			h_strich = h_strich-6.0
		c = (1- abs(2*(startL + i*(endL-startL)/steps)-1))*math.sqrt(startS + i*(endS-startS)/steps)
		x = c*(1- abs(h_strich%2-1))
		m = (startL + i*(endL-startL)/steps) -.5*c
		if (0.0 <= h_strich)and(h_strich<1.0):
			colormap.write('{'+str((c+m)*255.0)+','+str((x+m)*255.0)+','+str((0.0+m)*255.0)+'},\n')
		elif (1.0 <= h_strich)and(h_strich<2.0):
			colormap.write('{'+str((x+m)*255.0)+','+str((c+m)*255.0)+','+str((0.0+m)*255.0)+'},\n')
		elif (2.0 <= h_strich)and(h_strich<3.0):
			colormap.write('{'+str((0.0+m)*255.0)+','+str((c+m)*255.0)+','+str((x+m)*255.0)+'},\n')
		elif (3.0 <= h_strich)and(h_strich<4.0):
			colormap.write('{'+str((0.0+m)*255.0)+','+str((x+m)*255.0)+','+str((c+m)*255.0)+'},\n')
		elif (4.0 <= h_strich)and(h_strich<5.0):
			colormap.write('{'+str((x+m)*255.0)+','+str((0.0+m)*255.0)+','+str((c+m)*255.0)+'},\n')
		elif (5.0 <= h_strich)and(h_strich<=6.0):
			colormap.write('{'+str((c+m)*255.0)+','+str((0.0+m)*255.0)+','+str((x+m)*255.0)+'},\n')
	c = (1- abs(2*(startL)-1))*(startS)
	x = c*(1- abs(h_strich%2-1))
	m = (startL) -.5*c
	if (0.0 <= h_strich)and(h_strich<1.0):
		colormap.write('{'+str((c+m)*255.0)+','+str((x+m)*255.0)+','+str((0.0+m)*255.0)+'}};')
	elif (1.0 <= h_strich)and(h_strich<2.0):
		colormap.write('{'+str((x+m)*255.0)+','+str((c+m)*255.0)+','+str((0.0+m)*255.0)+'}};')
	elif (2.0 <= h_strich)and(h_strich<3.0):
		colormap.write('{'+str((0.0+m)*255.0)+','+str((c+m)*255.0)+','+str((x+m)*255.0)+'}};')
	elif (3.0 <= h_strich)and(h_strich<4.0):
		colormap.write('{'+str((0.0+m)*255.0)+','+str((x+m)*255.0)+','+str((c+m)*255.0)+'}};')
	elif (4.0 <= h_strich)and(h_strich<5.0):
		colormap.write('{'+str((x+m)*255.0)+','+str((0.0+m)*255.0)+','+str((c+m)*255.0)+'}};')
	elif (5.0 <= h_strich)and(h_strich<=6.0):
		colormap.write('{'+str((c+m)*255.0)+','+str((0.0+m)*255.0)+','+str((x+m)*255.0)+'}};')
	