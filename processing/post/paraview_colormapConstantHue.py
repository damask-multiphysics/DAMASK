#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

# This script is used to generate colormaps for paraview (www.paraview.org)
# The script writes 360 files. Each file contains one colormap.
# More information on the used colors space can be found at http://en.wikipedia.org/wiki/HSL_and_HSV
# written by M. Diehl, m.diehl@mpie.de

import math

print '******************************************************************************'
print '                   Write colormaps for paraview'
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
	colormap = open('colormap_' + str(h) + '.xml',"w")
	colormap.write('<ColorMap name = "' + str(h) + '" space = "RGB">\n')
	for i in range(0,steps+1):
		h_strich = h/60.0
		if(h_strich>6.0):
			h_strich = h_strich-6.0
		c = (1- abs(2*(startL + i*(endL-startL)/steps)-1))*math.sqrt(startS + i*(endS-startS)/steps)
		x = c*(1- abs(h_strich%2-1))
		m = (startL + i*(endL-startL)/steps) -.5*c
		if (0.0 <= h_strich)and(h_strich<1.0):
			colormap.write('<Point x="'+str(i)+'" o="1" r="'+str(c+m)+'" g="'+str(x+m)+'" b="'+str(0.0+m)+'"/>\n')
		elif (1.0 <= h_strich)and(h_strich<2.0):
			colormap.write('<Point x="'+str(i)+'" o="1" r="'+str(x+m)+'" g="'+str(c+m)+'" b="'+str(0.0+m)+'"/>\n')
		elif (2.0 <= h_strich)and(h_strich<3.0):
			colormap.write('<Point x="'+str(i)+'" o="1" r="'+str(0.0+m)+'" g="'+str(c+m)+'" b="'+str(x+m)+'"/>\n')
		elif (3.0 <= h_strich)and(h_strich<4.0):
			colormap.write('<Point x="'+str(i)+'" o="1" r="'+str(0.0+m)+'" g="'+str(x+m)+'" b="'+str(c+m)+'"/>\n')
		elif (4.0 <= h_strich)and(h_strich<5.0):
			colormap.write('<Point x="'+str(i)+'" o="1" r="'+str(x+m)+'" g="'+str(0.0+m)+'" b="'+str(c+m)+'"/>\n')
		elif (5.0 <= h_strich)and(h_strich<=6.0):
			colormap.write('<Point x="'+str(i)+'" o="1" r="'+str(c+m)+'" g="'+str(0.0+m)+'" b="'+str(x+m)+'"/>\n')       
	colormap.write('</ColorMap>')