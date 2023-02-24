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

import matplotlib.animation as animation
# from pylab import *
# import cmocean


""" *****************************************************************
"""

# dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref/"
# dirm = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref/"

dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_bvp"
dirm = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_bvp"


timeframe = "2h"
pdt = "/RIDGE_sco_1_%s_grid_T.nc" % (timeframe)
pdu = "/RIDGE_sco_1_%s_grid_U.nc" % (timeframe)
pdv = "/RIDGE_sco_1_%s_grid_V.nc" % (timeframe)
pdw = "/RIDGE_sco_1_%s_grid_W.nc" % (timeframe)
# pdw = "/RIDGE_ref_1_12h_grid_W.nc"
pmm = "/mesh_mask.nc"

save = 1 ; psave = "ridge" ; film = 1

########################################################
vmin = -2. ; vmax = 2. ; tskip=3
Ncolor=8
dv=(vmax-vmin)/Ncolor
"""
    Plot du test overflow - zps : z partial cells
    34 : time step
    101 : depth (z)
    3 : y dimension (1 cell embraced by two boundary layers)
    202 : x dimension
"""
dtu = nc4.Dataset(dir+pdu)
dt  = nc4.Dataset(dir+pdt)
mm = nc4.Dataset(dirm+pmm)

uu   = dtu.variables['u_vol'][::tskip,:,:,:]   # volumic transport
# toce = dt.variables['toce'] [::tskip,:,:,:]
nT,nK,nY,nX = np.shape(uu)
# nT=1

tmask = mm.variables['tmask'][0][:,:,:]
glamt = mm.variables['glamt'][0] ; gphit = mm.variables['gphit'][0]
glamu = mm.variables['glamu'][0]
umask = mm.variables['umask'][0] ; vmask = mm.variables['vmask'][0]
gdept = mm.variables['gdept_0'][0] ; gdepw = mm.variables['gdepw_0'][0]
e3w = mm.variables['e3w_0'][0]
mbathy= mm.variables['mbathy'][0,:,:]

print("domain size is (x,y) %dx%d with %d k-levels" % (nX,nY,nK))
midY = nY//2# np.where(np.abs(gphit[:,0])<=1E3)[0][0]
# attention en BVP le milieu est un V point
midX = np.where(np.abs(glamt[0,:])<=1E3)[0][0]

#######################################################
# streamfunction are UW points
yw = np.copy(gdepw[:,midY,:-1])
for i in range(nX-1):
    yw[:,i]=0.5*( gdepw[:,midY,i] + gdepw[:,midY,i+1] )
xu = yw*0.
for k in range(nK):
    xu[k,:]=glamu[midY,:-1]/1E3

# ... surrounded by T points
dz = np.copy(e3w[0,midY,:])
yt = np.zeros((nK+1,nX))
yt[0,:] = 0. ; yt[1:,:] = gdept[:,midY,:]
xt = yt*0. ; x1t = glamt[  midY,:]/1E3
for k in range(nK+1):
    xt[k,:]=x1t

########################################################
sfm = np.zeros((nT,nK,nX-1))
""" streamfunction """
for t in range(nT):
    sys.stdout.write(u"\u001b[1000D" + "processing psi [%3d/%3d]" % (t+1,nT))
    sys.stdout.flush()
    uum  = np.ma.masked_where(umask==0,uu[t])
    for i in range(nX-1):
        for k in range(nK):
            # k=0 surface - k=nK-1 last dot (masked)
            # there are as many u than psi
            sfm[t,k,i] += - np.sum(uum[k:,:,i])/1E6
sfm = np.ma.array(sfm)

# bigmask = 1-np.repeat(umask[np.newaxis,:,:,:],nT,axis=0)
# uu = np.ma.array(uu,mask=bigmask)
# for i in range(nX-1):
#     for k in range(nK):
#         sys.stdout.write(u"\u001b[1000D" + "processing psi [%3d/%3d]" % (1+i+k*(nX-1),nK*(nX-1)))
#         sys.stdout.flush()
#         sfm[:,k,i] += - np.sum(uu[:,k:,:,i],axis=(1,2))/1E6
# bigmask = 1-np.repeat(umask[np.newaxis,0,:,:],nT,axis=0)
# sfm = np.ma.array(sfm,mask=bigmask)


########################################################
# to have a look at the meshmask
# first bottom cell (masked - task[mbathy,:,:]-> 0)
zht = mbathy*0.
for ii in range(nX):
    for jj in range(nY):
        kk = mbathy[jj,ii]
        zht[jj,ii] = gdepw[kk,jj,ii]

########################################################


titlezer  = "Streamfunction"
palette = plt.get_cmap('RdBu_r',Ncolor)
optpcolor = {"vmin":vmin, "vmax":vmax, "cmap" : palette}
opthatch = {'facecolor':'grey', 'alpha' : 0.5, 'interpolate':True, 'step' : 'mid'}
levelsc = np.linspace(vmin,vmax, Ncolor+1)
optcontour = {'levels' : levelsc, 'colors':'black', 'linestyles' : "solid"}

fig, ax = plt.subplots(figsize=(10,6),dpi=200)

c  = ax.contour(xu,yw,sfm[0], **optcontour)
cf = ax.pcolormesh(xt, yt, sfm[0],**optpcolor)

cbar = plt.colorbar(cf)
cbar.set_label(r"Zonal streamfunction $\psi$ (Sv)")
# ax.fill_between(xu[0,:], strait(xw[0,:]), ridge(xw[0,:]),**opthatch)
# ... we would expect UW points to be the bottom of the basin
ax.fill_between(glamt[0,:]/1E3, zht[0,:], zht[midY,:], **opthatch)

if timeframe=="12h" :
    titlezer = '%02dd/%02dd \n'%((0.)*tskip/2.,nT*tskip/2.) # because 12h sortie
elif timeframe=="2h":
    titlezer = '%02dh/%02dh (%02dd/%02dd) \n'%((0.)*tskip*2.,nT*tskip*2.,(0.)*tskip/12.,nT*tskip/12.) # because 2h sortie
titlezer += "\nmin = %1.1f Sv   max = %1.1f Sv" % (np.nanmin(sfm[0]), np.nanmax(sfm[0]))
ax.set_title(titlezer, fontsize = 12, y = 1.02)

ax.set_ylim(5500,0)
ax.set_yticks([0,1000,2000,3000,4000,5000])
# ax.set_yticks([0,1000,2000,3000,4000,5000])
ax.set_ylabel("Z (m)")

ax.set_xlim(-1000,1000)
# ax.set_xlim(-2000,2000)
ax.set_xlabel("X (km)")

ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.yaxis.set_minor_locator(MultipleLocator(250))
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
plt.tight_layout()

def animate(i):
    """Set the data for the ith iteration of the animation."""
    global c,cf, sfm

    ax.collections = []
    c = ax.contour(xu,yw,sfm[i], **optcontour)
    cf = ax.pcolormesh(xt, yt, sfm[i], **optpcolor)
    #
    ax.fill_between(glamt[0,:]/1E3, zht[0,:], zht[midY,:], **opthatch)
    #
    if timeframe=="12h" :
        ptitle = '%02dd/%02dd \n'%((i)*tskip/2.,(nT-1)*tskip/2.) # because 12h sortie
    elif timeframe=="2h":
        ptitle = '%02dh/%02dh (%02dd/%02dd) \n'%((i)*tskip*2.,(nT-1)*tskip*2.,(i)*tskip/12.,(nT-1)*tskip/12.) # because 2h sortie
    ptitle += "min = %2.2f Sv   max = %2.2f Sv   " % (np.nanmin(sfm[i]), np.nanmax(sfm[i]))
    sys.stdout.write(u"\u001b[1000D" + "processing movie [%3d/%3d]" % (i+1,nT))
    sys.stdout.flush()
    ax.set_title(ptitle, fontsize = 12, y = 1.02)
    #
    return c, cf

if save:
    anim = animation.FuncAnimation(fig, animate, frames=nT, blit=False, repeat=False)
    writer = animation.writers['ffmpeg'](fps=4)
    anim.save('%s.mp4' % (psave), writer=writer, dpi=200)
    plt.close("all")
    print("\nsaving : %s" % psave)
else:
    anim = animation.FuncAnimation(fig, animate, frames=nT)
    plt.show()
