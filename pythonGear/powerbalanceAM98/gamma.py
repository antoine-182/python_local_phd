# -*- coding: utf-8 -*-

import netCDF4 as nc4
import numpy as np
import time, sys, os

import matplotlib.pyplot as plt
from pylab import *

""" *****************************************************************
"""
rn_abp = 0.01
D = 1.
gamma = D
l = 1.5*D

""" Tentative de polynome """
# def coeffM():
#     detA = gamma*l*(gamma-l)
#     za = 2*l*(gamma-1)/detA
#     zb = (l*l - gamma*l*l + l*l*l - l*gamma)/detA
#     zc = l
#     return((za,zb,zc))
# def invM():
#     a = np.array([ [gamma**2, gamma, 1.], \
#                    [    l**2,     l, 1.], \
#                    [     0. ,    0., 1.]  ])
#     ainv=np.linalg.inv(a)
#     print("%s" % np.allclose(np.dot(a, ainv), np.eye(3)) )
#     b = np.array([ [gamma - 1 ], \
#                  [        0.  ], \
#                  [          l ]  ])
#
#     return(np.dot(ainv, b))
# (a,b,c) = invM()
# def frpo(x,y):
#     return(y - a*x*x - b*x - c)

# """ Truc linéaire de polynome """
# def frpo(x,y):
#     return(y - l + x)

""" Cercle """
def frpo(x,y):
    # rajouter un max sur x/y
    return( (y-l)**2 + (x-l)**2 - l**2)

xlist = np.arange(0.,4*D,D/20.)
ylist = np.copy(xlist)
X,Y = np.meshgrid(xlist,xlist)

rpo = frpo(X,Y)
one=np.ones(np.shape(rpo))
mrpo=one
for i in range(len(xlist)):
    for j in range(len(ylist)):
        z1d = 1 - rpo[i,j]
        mrpo[i,j] =np.nanmin( (np.nanmax((rn_abp,z1d)),1.) )


palette = plt.get_cmap('binary')
fig, ax = plt.subplots(dpi=200)
im=ax.pcolormesh(X,Y,mrpo, alpha =0.9, cmap=palette)
plt.colorbar(im)
ax.axvline(x= 0, linestyle = '-', color = 'red', lw = 1.)
ax.axhline(y= 0, linestyle = '-', color = 'red', lw = 1.)
# plt.plot(xlist,a*xlist*xlist+b*xlist+c,"red")
# plt.plot(xlist,l-xlist,"red")
plt.show()
