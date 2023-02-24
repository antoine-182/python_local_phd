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

Listdt_up = [[0, "", "/fluxform"],
             [-1, "", "/fluxform"],
             [-1, "", "/vectformPS"],
             [-1, "", "/vectformNS"]]
Listdt_do = [[0, "", "x"],
             [-1, "", "/freeslip"],
             [-1, "", "/divNS"],
             [-1, "", "/fauxFS"]]

Listdt = [Listdt_up,Listdt_do]

pdF = "/CORNER_1d_00010101_00010130_grid_F.nc"
# pdw = "/RIDGE_ref_1_12h_grid_W.nc"
pmm = "/mesh_mask.nc"

vmax = 1e-5 ; tskip=1
Ncolor=7
dv=2*vmax/Ncolor
dx = 5. # km

########################################################

palette = plt.get_cmap('RdBu_r',Ncolor)
cticks = np.arange(-vmax,vmax+dv,dv)
optpcolor = {"vmin":-vmax, "vmax":vmax, "cmap" : palette, "alpha" : 0.9}
optcontour = {"vmin":-vmax,"vmax":vmax, "levels":cticks, "linewidths" :0.1, "colors":('k',),"linestyles" : "solid"}


save = 1 ; psave = "figure_corner.png" ; dpi = 200
NI = 2 ; NJ = 4 # nb of rows (NI) and cols (NJ)
fig, ax = plt.subplots(NI,NJ, figsize=(NJ*2,NI*2), dpi = dpi, squeeze=False)
# settled per rows (i=0) j=0 j=1 j=2
#                  (i=1) j=0 j=1 j=2
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

""" *************************** DATA
"""
########################################################
#init
N=len(Listdt)
data=np.zeros((nY,nX))

for i in range(NI):
    for j in range(NJ):
        frame, tit,dirname = Listdt[i][j]
        print("[%d,%d] : %s" % (i,j,dirname))
        if dirname == 'x':
            ax[i][j].axis('off')
        else :
            dtf = nc4.Dataset(dir+dirname+pdF)
            data[:,:] = dtf.variables['relvor'][frame,:,:]   # zeta
            data = np.ma.masked_where(fmask_vorlatF==0,data)
            # pcolormesh
            cf = ax[i][j].pcolormesh(glamt, gphit,data ,**optpcolor)
            ax[i][j].set_title(tit, y=1.04, fontsize = 8)
            # contour
            ct= ax[i][j].contour(glamf, gphif, data, **optcontour)
            # ax[i][j].clabel(ct, fmt='%e', colors='k', fontsize=8)


""" *************************** PRETTY
"""
for i in range(NI):
    for j in range(NJ):
        vorlat, tit,dirname = Listdt[i][j]
        if dirname == 'x':
            continue
        # print(tickszer)
        ax[i][j].set_xlim(-dx/2,500+dx/2)
        ax[i][j].set_ylim(-dx/2,500+dx/2)
        ax[i][j].set_xticklabels([]) ; ax[i][j].set_yticklabels([])
        #
        if (i==NI-1):
            ax[i][j].set_xticks([0,200,400])
            ax[i][j].set_xticklabels(["0","200","400"])
            ax[i][j].set_xlabel("X (km)")
        #
        if (j==NJ-1):
            ax[i][j].set_yticks([0,200,400])
            ax[i][j].set_yticklabels(["0","200","400"])
            ax[i][j].set_ylabel("Y (km)")
            ax[i][j].yaxis.set_label_position("right")
        #
        ax[i][j].patch.set_color('0.8')
        ax[i][j].xaxis.set_minor_locator(MultipleLocator(50))
        ax[i][j].yaxis.set_minor_locator(MultipleLocator(50))
        ax[i][j].tick_params(axis = "y", which = 'both', width=1.2, labelsize = 8, pad = 3, left = True, right = True)
        ax[i][j].tick_params(axis = 'x', which = 'both', width=1.2, labelsize = 8, pad = 3, bottom = True, top = True)
        ax[i][j].tick_params(which='minor',length = 3)
        ax[i][j].tick_params(which='major',length = 4)
        ax[i][j].set_aspect(aspect='equal') # data coordinate 'equal'

#fig.subplots_adjust(hspace = 0., wspace = 0.2, right=0.85)
#plt.subplots_adjust(hspace = )
plt.subplots_adjust(wspace=-0.5, hspace=0.2)
#fig.tight_layout(pad=.2)

cticks = np.array([-vmax,-vmax+dv*2,vmax-dv*2,vmax])

cbar=plt.colorbar(cf, ax=ax[:], shrink=0.5, location = "bottom",
                  aspect = 40, fraction=0.05,  extend = 'both', pad = 0.2) #
cbar.set_ticks(cticks)
cbar.set_ticklabels(["%1.e" % s for s in cticks])
cbar.set_label(r"Relative vorticity $\zeta/h \ [\mathrm{m}^{-1}\mathrm{s}^{-1}]$")

""" *************************** SAVING
"""

if save :
    print("\nsaving : %s" % psave)
    fig.savefig(psave, dpi = dpi)
    plt.close()
    print("figure closed")
plt.show()
