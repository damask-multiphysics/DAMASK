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
print 'Suitable for datasets running from negative to positive values.'
print 'The colors are described using the HSL model.'
print 'Each of the 360 generated colormaps uses two values of (H)ue.'
print 'The colormaps runs at constant H from given (L)ightness and (S)aturation'
print 'to given L_min and S_min and goes to H+180° (with given L and S)'
print 'Suitable values: L = L_min =.5, S = 1, and S_min=0,'
print '******************************************************************************'
print ''
startL = float(raw_input('Please enter start value for (L)ightness: '))
endL = float(raw_input('Please enter minimum value for L: '))
startS = float(raw_input('Please enter start value for (S)aturation: '))
endS = float(raw_input('Please enter minimum value for S: '))
steps = int(raw_input('Please enter steps/resolution: '))
steps = steps/2
for h in range(0,360):
	colormap = open('colormap_' + str(h).zfill(3) + '.xml',"w")
	colormap.write('<ColorMap name = "Complementary Color ' + str(h).zfill(3) + '" space = "RGB">\n')
	i=0
	h_strich = h/60.0
	if(h_strich>6.0):
			h_strich = h_strich-6.0
	for j in range(0,steps+1):
		# let L run linearly from 0 to 1, let S be the square root of a linear list from 0 to 1
		c = (1- abs(2*(startL - j*(startL-endL)/steps)-1))*(startS - j*(startS-endS)/steps)
		x = c*(1- abs(h_strich%2-1))
		m = (startL - j*(startL-endL)/steps) -.5*c		
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
		i = i+1	
	h_strich = (h+180.0)/60.0
	if(h_strich>6.0):
		h_strich = h_strich-6.0
	for j in range(1,steps+1):
		c = (1- abs(2*(endL+j*(startL-endL)/steps)-1))*(endS+j*(startS-endS)/steps)
		x = c*(1- abs(h_strich%2-1))
		m = (endL+j*(startL-endL)/steps) -.5*c
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
		i=i+1			
	colormap.write('</ColorMap>')