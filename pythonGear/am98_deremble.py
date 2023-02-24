# -*- coding: utf-8 -*-

import netCDF4 as nc4
import numpy as np
import time, sys, os

# one way to do
# sys.path.append(os.path.abspath('vacumm-3.4.0/'))

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, FixedLocator, FixedFormatter,
                               NullLocator)
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as animation
# from pylab import *
# import cmocean


""" *****************************************************************
"""

dir = "/Users/gm/Documents/nemo/dev_r12527_Gurvan_ShallowWater/cfgs/CORNER/EXP_article"

Listdt_up = [[0, "", "/init"],
             [1, "", "/fluxform"],
             [0, "", "/vectformPS"],
             [0, "", "/vectformNS"]]
Listdt_do = [[0, "", "x"],
             [0, "", "/fauxFS"],
             [1, "", "/divNS"],
             [0, "", "/fauxFS"]]
Listdt = [Listdt_up,Listdt_do]

pdF = "/CORNER_1d_00010101_00010230_grid_F.nc"
# pdw = "/RIDGE_ref_1_12h_grid_W.nc"
pmm = "/mesh_mask.nc"

save = 0 ; psave = "deremble" ; film = 1

vmax = 1e-5 ; tskip=1
Ncolor=7
dv=2*vmax/Ncolor
dx = 5. # km
########################################################

mm = nc4.Dataset(dir+pmm)
tmask = mm.variables['tmask'][0,0]
nY,nX = np.shape(tmask)
fmask_vorlatT = np.zeros((nY,nX)) ; fmask_vorlatF = np.zeros((nY,nX))
for i in range(nX-1):
    for j in range(nY-1):
        fmask_vorlatT[j,i] = tmask[j+1,i] * tmask[j+1,i+1] \
                           * tmask[j  ,i] * tmask[j  ,i+1]
        fmask_vorlatF[j,i] = tmask[j+1,i] * tmask[j+1,i+1] \
                           * tmask[j  ,i] * tmask[j  ,i+1]
        # exception vorlat
        if ((tmask[j+1,i] + tmask[j+1,i+1] \
            + tmask[j  ,i] + tmask[j  ,i+1])) == 3 :
            fmask_vorlatT[j,i] = 1.

glamt = mm.variables['glamt'][0]/1E3 ; gphit = mm.variables['gphit'][0]/1E3 # corner of cells
glamf = mm.variables['glamf'][0]/1E3 ; gphif = mm.variables['gphif'][0]/1E3 # F points

########################################################

palette = plt.get_cmap('RdBu_r',Ncolor)
cticks = np.arange(-vmax,vmax+dv,dv)
optpcolor = {"vmin":-vmax, "vmax":vmax, "cmap" : palette, "alpha" : 0.9}
optcontour = {"vmin":-vmax,"vmax":vmax, "levels":cticks, "linewidths" :0.1, "colors":('k',),"linestyles" : "solid"}

# for _ in range(len(cticks)):
#     if _%2==1:
#         cticks[_]=np.nan

fig, ax = plt.subplots(1,4, figsize=(3*4,1*4),
                       # constrained_layout=True,
                       dpi=200)

########################################################
#init
N=len(Listdt)
nT=np.int64(60/tskip)
zeta=np.zeros((N,nT,nY,nX))
for n in range(N):
    vorlat, tit,dirname = Listdt[n]
    dtf = nc4.Dataset(dir+dirname+pdF)
    zeta[n,:,:,:] = dtf.variables['relvor'][::tskip,:,:]   # zeta
    if vorlat :
        data = np.ma.masked_where(fmask_vorlatT==0,zeta[n,0])
    else :
        data = np.ma.masked_where(fmask_vorlatF==0,zeta[n,0])
    # pcolormesh
    cf = ax[n].pcolormesh(glamt, gphit,data ,**optpcolor)
    ax[n].set_title(tit, y=1.04, fontsize = 8)
    # contour
    ct= ax[n].contour(glamf, gphif, data, **optcontour)
    # ax[n].clabel(ct, fmt='%e', colors='k', fontsize=8)

# Cosmetic
for n in range(N):
    vorlat, tit,dirname = Listdt[n]
    # print(tickszer)
    ax[n].set_xlim(-dx/2,500+dx/2)
    ax[n].set_ylim(-dx/2,500+dx/2)
    #
    ax[n].set_xlabel("X (km)", fontsize=8)
    #
    if (n==0):
        ax[n].set_ylabel("Y (km)", fontsize=8)
    if (n>0):
        ax[n].set_yticklabels([])
    #
    ax[n].patch.set_color('0.8')
    ax[n].xaxis.set_minor_locator(MultipleLocator(50))
    ax[n].yaxis.set_minor_locator(MultipleLocator(50))
    ax[n].tick_params(axis = "y", which = 'both', width=1.2, labelsize = 8, pad = 3, left = True, right = True)
    ax[n].tick_params(axis = 'x', which = 'both', width=1.2, labelsize = 8, pad = 3, bottom = True, top = True)
    ax[n].tick_params(which='minor',length = 3)
    ax[n].tick_params(which='major',length = 4)
    ax[n].set_aspect(aspect='equal') # data coordinate 'equal'

plt.subplots_adjust(hspace=0.8)

cbar=plt.colorbar(cf, ax=ax[:], shrink=0.5, location = "bottom",
                  aspect = 40, fraction=0.05,  extend = 'both', pad = 0.15) #
cbar.set_ticks(cticks)
cbar.set_ticklabels(["%1.e" % s for s in cticks])
cbar.set_label(r"Relative vorticity $\zeta/h \ [\mathrm{m}^{-1}\mathrm{s}^{-1}]$")



def animate(i):
    """Set the data for the ith iteration of the animation."""
    for nn in range(N) :
        vorlat,_,_ = Listdt[nn]
        if vorlat :
            data = np.ma.masked_where(fmask_vorlatT==0,zeta[nn,i])
        else :
            data = np.ma.masked_where(fmask_vorlatF==0,zeta[nn,i])
        ax[nn].collections = []
        cf = ax[nn].pcolormesh(glamt, gphit, data, **optpcolor)
        ct = ax[nn].contour(glamf, gphif, data, **optcontour)
        # ax[nn].clabel(ct, fmt='%e', colors='k', fontsize=8)
        #
        sys.stdout.write(u"\u001b[1000D" + "processing movie %d/%d[%3d/%3d]" % (nn,N,i+1,nT))
        sys.stdout.flush()

    return cf, ct

# ptitle = '%02dd/%02dd'%((i+1)*tskip/2.,nT*tskip/2.) # because 12h sortie
# ax[nn].set_title(ptitle, fontsize = 12, y = 1.02)

if save:
    anim = animation.FuncAnimation(fig, animate, frames=nT, blit=False, repeat=False)
    writer = animation.writers['ffmpeg'](fps=16)
    anim.save('%s.mp4' % (psave), writer=writer, dpi=200)
    plt.close("all")
    print("\nsaving : %s" % psave)
else:
    anim = animation.FuncAnimation(fig, animate, frames=nT)
    plt.show()
