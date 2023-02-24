# -*- coding: utf-8 -*-

import time, sys, os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy import interpolate

""" Teste un raccordement C1
    Fonctionne que pour l1=l2
"""

""" slope """
beta = 20. # m/km
alpha = 500. # m
l1 = 5 ; l2 = 5

x = np.arange(0,200,0.5) ; x0 = 30 ; x1 = x0 + 1500./beta # km
ix0 = np.argwhere(x==x0)[0][0] ; ix1 = np.argwhere(x==x1)[0][0]
zht = 2000. + x*0. ; zht[:ix0] = 500. ; zht[ix0:ix1]= beta*(x[ix0:ix1]-x0) + alpha

""" spline x0 """
za1 = (l1-l2)/(l1+l2)
za2 = 1/(l1*l1 + l2*l2 - l1*l2)
za3 = beta/4.
a = za3 * za2 * za1

zb1 = beta/(2*(l1+l2))
zb2 = 3*a*(l1-l2)
b = zb1 + zb2

zc1 = beta*l1/(l1+l2)
zc2 = 3*a*l1*l2
c = zc1 - zc2

# zd1 = 2*a*(l1*l1*l2*l2)/(l1-l2)
zd1 = 2*(l1*l1*l2*l2)/(l1+l2)
d = alpha + zd1*za2*za3

fl = np.poly1d(np.array(( a,b,c,d )) )
polx = fl(x-x0)

""" spline x1 """
beta = 20. # m/km
alpha = 2000. # m
l1 = -l1 ; l2 = -l2

za1 = (l1-l2)/(l1+l2)
za2 = 1/(l1*l1 + l2*l2 - l1*l2)
za3 = beta/4.
a = za3 * za2 * za1

zb1 = beta/(2*(l1+l2))
zb2 = 3*a*(l1-l2)
b = zb1 + zb2

zc1 = beta*l1/(l1+l2)
zc2 = 3*a*l1*l2
c = zc1 - zc2

# zd1 = 2*a*(l1*l1*l2*l2)/(l1-l2)
zd1 = 2*(l1*l1*l2*l2)/(l1+l2)
d = alpha + zd1*za2*za3

fl = np.poly1d(np.array(( a,b,c,d )) )
poly = fl(x-x1)

""" plot """
plt.plot(x,zht,  '-', c = "black", marker = '.')
plt.plot(x,polx, '--', c = "red", marker = 'o')
plt.plot(x,poly, '--', c = "royalblue", marker = 'o')
plt.vlines(x0-np.abs(l1), 0, 2000)
plt.vlines(x0+np.abs(l2), 0, 2000)
plt.vlines(x1-np.abs(l2), 0, 2000)
plt.vlines(x1+np.abs(l1), 0, 2000)
plt.xlabel("lenght (km)")
plt.ylabel("depth (m)")
plt.xlim(1,200)
plt.ylim(2200.,200.)
plt.show()
