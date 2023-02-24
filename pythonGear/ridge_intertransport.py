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

pdu = "/RIDGE_sco_1_12h_grid_U.nc"

dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref"
dirm = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref"

Listdt = [[1 ,"1° U strait","/ref1U" ],
          [1 ,"1° V strait","/ref1F" ]]


save = 0 ; psave = "transport.png"
utype = 2  # U>0
           # psi
           # U>0 + U<0

N = len(Listdt)
if (utype == 0 or utype == 1) :
    Q = np.zeros((720, N))
elif utype == 2 :
    Q = np.zeros((720,3, N))

rn_h_ridge = 2800.

for n in range(N):
    nn,lin,name = Listdt[n]
    dtu = nc4.Dataset(dir +name + pdu)
    #
    pmm = "/mesh_mask.nc"
    mm  = nc4.Dataset(dir +name + pmm)
    #
    print("data dataframe %d/%d" % (n+1,N))
    tu  = dtu.variables['u_vol'] # transport penalised
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
    zone = (a==1)*(b>rn_h_ridge)  # sous le ridge
    # zone = (a==1)*(b<rn_h_ridge)  # au dessus du ridge
    #
    nTok = np.min((720,nT))
    if (utype == 0) :
        for t in range(nTok):
            ### intégration simple
            #uum = um[t][(a==1)*(b>rn_h_ridge)]
            #Q[t,n] = np.sum(uum)
            ### U>0
            uum = um[t][zone]
            Q[t,n] = np.sum(uum[uum>0])
    elif (utype == 2) :
        for t in range(nTok):
            ### intégration simple
            uum = um[t][zone]
            Q[t,0,n] = np.sum(uum[uum>0])
            Q[t,1,n] = np.sum(uum[uum<0])
            Q[t,2,n] = np.sum(uum)
    elif (utype==1):
        for t in range(nTok):
            ### fonciton de courant
            uum = um[t][::-1,:]
            psi=depw*0
            for k in range(1,len(depw)) :
                z = depw[k]
                psi[k] = np.sum(uum[zdep>z])
            Q[t,n] = np.max(psi)

## intégration simple
##uum = um[t][(a==1)*(b>rn_h_ridge)]
##Q[t,n] = np.sum(uum)
## amélioration possible
#uum = um[t][::-1,:] ;
#psi=depw*0
#for k in range(1,len(depw)) :
#    z = depw[k]
#    psi[k] = np.sum(uum[zdep>z])
#Q[t,n] = np.max(psi)

Q/=1e6


# In[ ]:


""" figure """
# fig, ax = plt.subplots(figsize=(8, 6),dpi=200)
fig, ax = plt.subplots()
timelist = (np.arange(0,nT,1)+1)/2
# at each front
for n in range(N):
    nn,lin,name = Listdt[n]
    print("plot dataframe %d/%d" % (n+1,N))
    dtu = nc4.Dataset(dir +name + pdu)
    #
    tu  = dtu.variables['u_vol'] # transport penalised
    nT,_,_,_ = np.shape(tu)
    nTok = np.min((720,nT))
    #
    timelist = (np.arange(0,nTok,1)+1)/2
    if (utype == 0 or utype == 1) :
        line, = ax.plot(timelist, Q[:nTok,n],
                        label=lin,
                        alpha=1., linewidth=1.5)
    elif utype == 2 :
        line, = ax.plot(timelist, Q[:nTok,2,n],
                        label=lin,
                        alpha=1., linewidth=1.5)
        ax.plot(timelist, Q[:nTok,0,n],   # >0
                label="_no_legend_",
                color = line.get_color(), linestyle = "--",
                alpha=1., linewidth=1.5)
        ax.plot(timelist, -Q[:nTok,1,n],   # <0
                label="_no_legend_",
                color = line.get_color(), linestyle = ":",
                alpha=1., linewidth=1.5)


ax.grid(which='major', linestyle='-', linewidth='0.3', color='black')
ax.grid(which='minor', linestyle=':', linewidth='0.3', color='black')
# ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_major_locator(MultipleLocator(60.))
ax.xaxis.set_minor_locator(MultipleLocator(30))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
ax.set_xlim(0,180)
ax.set_xlabel("Time in days")
# ax.set_yscale('log')
# ax.set_ylim(1e-2,2e1)
ax.set_ylabel("Volume Transport (Sv)")
ax.legend()

fig.add_subplot(ax)

titlezer = "Transport across the ridge\n"
# titlezer += "min = %2.2f kg/m3   max = %2.2f kg/m3   " % (np.min(toce[toce!=0]), np.max(toce[toce!=0]))
ax.set_title(titlezer, fontsize = 12, y = 1.02)
plt.tight_layout(rect=[0,0,1,0.95])


if save :
    print("saving : %s" % psave)
    fig.savefig(psave)
    plt.close()
else :
    plt.show()
