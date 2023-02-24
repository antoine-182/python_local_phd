
# -*- coding: utf-8 -*-

import netCDF4 as nc4
import numpy as np
import time, sys

# one way to do
# sys.path.append(os.path.abspath('vacumm-3.4.0/'))

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, FixedLocator, FixedFormatter,
                               NullLocator)
import matplotlib.gridspec as gridspec

import matplotlib.animation as animation
# from pylab import *
# import cmocean
import types

dx = 1000. 
gridx = np.arange(0,21*dx,dx) 
Ni = len(gridx) 

def shapiro_1d(tab,n=1):
    z2d = np.copy(tab)
    for _ in range(n) :
        for i in range (1,len(tab)-1) :
            z2d[i] = 0.25 * tab[i-1] + 0.5 * tab[i] + 0.25 * tab[i+1]
        tab=z2d
    return(tab)


rn_abp = 1e-3 ; smo = 0 ; minmoy = 1
mid = Ni//2
a = np.ones((Ni))
a[:mid] = rn_abp

rpot = shapiro_1d(a,smo)

rpou=np.copy(rpot)
if minmoy==1 : 
    for i in range(1,Ni-1):
        rpou[i] = 0.5*(rpot[i]+rpot[i+1])
else : 
    for i in range(1,Ni-1):
        rpou[i] = np.min(rpot[i],rpot[i+1])

plt.plot(rpot,label=r"$\phi_t$")
plt.plot(rpou,label=r"$\phi_u$")
plt.plot(rpot/rpou,label=r"$\phi_u/\phi_t$")
plt.legend()
plt.show()
