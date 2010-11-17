#!/usr/bin/python

#
#    Copyright 2010 Andrew Walker
#
#    This file is part of the Dislocator package.
#
#    The Dislocator package is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    The Dislocator package is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with the Dislocator package.  If not, see <http://www.gnu.org/licenses/>.
#


import matplotlib
matplotlib.use('PS')
from pylab import *
from numpy import *
import Image
import sys


im = Image.open('./001_background.jpg')
dpi = rcParams['figure.dpi']
figsize = im.size[0]/dpi, im.size[1]/dpi

min_x = -10.0
min_y = -12.8
max_x = 12.7
max_y = 12.6

mv_x = 0.1 # 0.6
mv_y = 0.1 #0.55

x_pos = []
y_pos = []
x_len = []
y_len = []

# build list of all URLs found in standard input
for line in sys.stdin:
        if line.startswith('ARROW'):
                #print line.rstrip()
                arrow = line.split()
		if ( (float(arrow[1]) < max_x) and (float(arrow[2]) < max_y) and (float(arrow[1]) > min_x) and (float(arrow[2]) > min_y) and (float(arrow[3]) < 2.0) and (float(arrow[3]) > -2.0) ):
          		x_pos.append((float(arrow[1])-min_x+mv_x)*0.25*dpi)
               		y_pos.append((float(arrow[2])-min_y+mv_y)*0.25*dpi)
               		x_len.append(float(arrow[5])*0.19*dpi)
               		y_len.append(float(arrow[6])*0.19*dpi)


#1
figure(figsize=figsize)
bgim = imshow(im, origin='lower')
axis('off')



Q = quiver(x_pos, y_pos, x_len, y_len, pivot='middle', units='x', scale=1)

#l,r,b,t = axis()
#dx, dy = r-l, t-b
#axis([l-0.05*dx, r+0.05*dx, b-0.05*dy, t+0.05*dy])

#show()
savefig('001_vit.ps')
