#!/usr/bin/env python
# coding: utf-8

# In[1]:


# -*- coding: utf-8 -*-

import netCDF4 as nc4
import numpy as np
import time, sys, os

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, FixedLocator, FixedFormatter,
                               NullLocator)
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
# from pylab import *
# import cmocean

#####################################################



dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref/sanity3"
dirm = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref/sanity3"

pmm = "/mesh_mask.nc"
pdu = "/RIDGE_sco_1_2h_grid_U.nc" ; to_day = 12.

save = 1 ; psave = "transport.png"

N = 2 # above and below
Q = np.zeros((720,3, N))

rn_h_ridge = 2800. ;  rn_hdepth = 5500.

dtu = nc4.Dataset(dir + pdu)
#
mm  = nc4.Dataset(dir + pmm)
#
tu  = dtu.variables['u_vol'] # transport (penalised)
nT,nK,nY,nX = np.shape(tu)
glamt = mm.variables['glamt'][0] ; gphit = mm.variables['gphit'][0]
umask = mm.variables['umask'][0]
gdept = mm.variables['gdept_0'][0] ; gdepw1d = mm.variables['gdepw_1d'][0] ;
#
dy = gphit[1,0] - gphit[0,0]
midY = np.where(np.abs(gphit[:,0])<=dy)[0][0]
midX = np.where(np.abs(glamt[0,:])<=1E3)[0][0]
# strait
um = tu[:,:,:,midX] ; depw = gdepw1d[::-1]
zdep = gdept[::-1,:,midX]   # profondeur
a = umask[:,:,midX] ; b = gdept[:,:,midX]

#
for n in range(N):
    if n == 0:
        zone = (a==1)*(b>rn_h_ridge)  # sous le ridge
    elif n==1:
        zone = (a==1)*(b<rn_h_ridge)  # au dessus du ridge
    for t in range(nT):
        ### intégration simple
        uum = um[t][zone]
        Q[t,0,n] = np.sum(uum[uum>0])
        Q[t,1,n] = np.sum(uum[uum<0])
        Q[t,2,n] = np.sum(uum)
Q/=1e6

""" figure """
# fig, ax = plt.subplots(figsize=(8, 6),dpi=200)
fig, ax = plt.subplots(figsize=(8,4),dpi=200)
timelist = (np.arange(0,nT,1)+1)/to_day
# at each front

for n in range(N):
    if n==0:
        namelabel = "below the ridge (%dm-%dm)" % (rn_h_ridge,rn_hdepth)
    elif n==1:
        namelabel = "above the ridge (0-%dm)" % (rn_h_ridge)
    line, = ax.plot(timelist, Q[:,2,n],
                    label=namelabel,
                    alpha=1., linewidth=1.5)
    ax.plot(timelist, Q[:,0,n],   # >0
            label="_no_legend_",
            color = line.get_color(), linestyle = "--",
            alpha=1., linewidth=1.5)
    ax.plot(timelist, Q[:,1,n],   # <0
            label="_no_legend_",
            color = line.get_color(), linestyle = ":",
            alpha=1., linewidth=1.5)


ax.grid(which='major', linestyle='-', linewidth='0.3', color='black')
ax.grid(which='minor', linestyle=':', linewidth='0.3', color='black')
# ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_major_locator(MultipleLocator(15.))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
# ax.set_xlim(0,180)
ax.set_xlabel("Days")
# ax.set_yscale('log')
# ax.set_ylim(1e-2,2e1)
ax.set_ylabel("Volume Transport (Sv)")
ax.legend()

fig.add_subplot(ax)

titlezer = "Zonal Transport"
# titlezer += "min = %2.2f kg/m3   max = %2.2f kg/m3   " % (np.min(toce[toce!=0]), np.max(toce[toce!=0]))
ax.set_title(titlezer, fontsize = 12, y = 1.02)
plt.tight_layout(rect=[0,0,1,0.95])


if save :
    print("saving : %s" % psave)
    fig.savefig(psave)
    plt.close()
else :
    plt.show()
