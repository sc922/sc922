# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 11:08:08 2017

@author: sc922
"""

from numpy import genfromtxt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.interpolate import interp1d

plt.close('all')

CIE = genfromtxt('ciexyz31.csv', delimiter=',')
wl_wave = genfromtxt('wl_wv.txt')
Spectra = genfromtxt('ANToM_40C.txt')

"""
CIE_Sensitivity = plt.figure()
plt.plot(CIE[:,0],CIE[:,1])
plt.plot(CIE[:,0],CIE[:,2])
plt.plot(CIE[:,0],CIE[:,3])
plt.ylabel('Intensity')
plt.xlabel('wavelength [nm]')
plt.show()
"""

Spectra_Fig = plt.figure()
plt.plot(wl_wave[:-1],Spectra)
plt.ylim([0,0.012])
plt.xlim([420,800])
plt.ylabel('Intensity')
plt.xlabel('wavelength [nm]')
plt.show()

CIE_red = interp1d(CIE[:,0],CIE[:,1])
CIE_grn = interp1d(CIE[:,0],CIE[:,2])
CIE_blu = interp1d(CIE[:,0],CIE[:,3])

Red_value = sum(CIE_red(wl_wave[20:529])*Spectra[20:529])
Grn_value = sum(CIE_grn(wl_wave[20:529])*Spectra[20:529])
Blu_value = sum(CIE_blu(wl_wave[20:529])*Spectra[20:529])

print Red_value, Grn_value, Blu_value

rR =  (3.2406*Red_value 	-1.5372*Grn_value 	-0.4986*Blu_value )
gG =  (-0.9689*Red_value 	+1.8758*Grn_value 	+0.0415*Blu_value )
bB =  (0.05575*Red_value 	-0.20404*Grn_value 	+1.05731*Blu_value)
	
if (rR > 0.0031308):
   rR = 1.055 * ( rR ** ( 1 / 2.4 ) ) - 0.055
else:                     
   rR = 12.92 * rR
				
if ( gG  > 0.0031308 ): 
   gG  = 1.055 * ( gG  ** ( 1 / 2.4 ) ) - 0.055
else:                   
   gG = 12.92 * gG
		
if ( bB > 0.0031308 ):
   bB = 1.055 * ( bB ** ( 1 / 2.4 ) ) - 0.055
else:                     
   bB = 12.92 * bB

if (rR > 1):
   r = 1
elif (rR <0):
   r = 0
else:
   r = rR

if (gG > 1):
   g = 1
elif (gG  <0):
   g = 0
else:
   g = gG 	

if (bB > 1):
   b = 1
elif (bB  <0):
   b = 0
else:
   b = bB
   
print r*255,g*255,b*255

class Palette:

    def __init__(self):
        self.palette = []

    def __call__(self, r, g, b):
        # map rgb tuple to colour index
        rgb = r, g, b
        try:
            return self.palette.index(rgb)
        except:
            i = len(self.palette)
            if i >= 256:
                raise RuntimeError, "all palette entries are used"
            self.palette.append(rgb)
            return i

    def getpalette(self):
        # return flattened palette
        palette = []
        for r, g, b in self.palette:
            palette = palette + [r, g, b]
        return palette
        
rgb = Palette()

#plt.show()

ax2 = Spectra_Fig.add_subplot(111)
ax2 = Spectra_Fig.add_axes([0.6,0.6,0.25,0.25])
ax2.xaxis.set_visible(False)
ax2.yaxis.set_visible(False)
ax2.set_zorder(1000)
ax2.patch.set_alpha(1)
ax2.patch.set_color([r, g, b])
"""
#"savefig('rect1.png', dpi=90, bbox_inches='tight')

im = Image.new("P", (400, 400), rgb(0, 0, 0))

d = ImageDraw.ImageDraw(im)
d.setfill(1)

d.setink(rgb(255, 0, 0))
d.polygon((0, 0, 0, 400, 400, 400))

d.setink(rgb(255, 153, 0))
d.rectangle((100, 100, 300, 300))

d.setink(rgb(255, 255, 0))
d.ellipse((120, 120, 280, 280))

im.putpalette(rgb.getpalette())

im.save("out.gif")


from graphics import *
win = GraphWin()
rect = Rectangle(Point(20, 10), pt)
rect.draw(win)
"""