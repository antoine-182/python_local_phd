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

# dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref/1ubs"
# dirm = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref/1ubs"

dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_bvp"
dirm = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_bvp"


timeframe = "2h"
pdt = "/RIDGE_sco_1_%s_grid_T.nc" % (timeframe)
pdu = "/RIDGE_sco_1_%s_grid_U.nc" % (timeframe)
pdv = "/RIDGE_sco_1_%s_grid_V.nc" % (timeframe)
pdw = "/RIDGE_sco_1_%s_grid_W.nc" % (timeframe)
# pdw = "/RIDGE_ref_1_12h_grid_W.nc"
pmm = "/mesh_mask.nc"
save = 1 ; psave = "ridge" ; film = 1 ; tskip=9

########################################################

dtu = nc4.Dataset(dir+pdu)
dt  = nc4.Dataset(dir+pdt)
mm = nc4.Dataset(dir+pmm)

tmask = mm.variables['tmask'][0][:,:,:]
nK,nY,nX = np.shape(tmask)

glamt = mm.variables['glamt'][0] ; gphit = mm.variables['gphit'][0]
glamu = mm.variables['glamu'][0]  ; gphiv = mm.variables['gphiv'][0]
umask = mm.variables['umask'][0] ; vmask = mm.variables['vmask'][0]
gdept = mm.variables['gdept_0'][0] ; gdepw = mm.variables['gdepw_0'][0]
e3w = mm.variables['e3w_0'][0]
mbathy= mm.variables['mbathy'][0,:,:]

uu   = dtu.variables['uoce'][::tskip,:,:,:]   # cross speed
toce = dt.variables['toce'][::tskip,:,:,:]
nT,_,_,_ = np.shape(toce)
# nT=1

print("domain size is (x,y) %dx%d with %d k-levels" % (nX,nY,nK))
midY = nY//2  # np.where(np.abs(gphit[:,0])<=1E3)[0][0]
# attention en BVP le milieu est un V point
midX = np.where(np.abs(glamt[0,:])<=1E3)[0][0]

#######################################################
# density points are T,T points
yt = np.zeros((nK,nY))
yt[:,:] = gdept[:,:,midX]
xt = yt*0. ; x1t = gphit[:,midX]/1E3
for k in range(nK):
    xt[k,:]=x1t

# Speed points are T,T points
# ... surrounded by UW
yuw = np.zeros((nK+1,nY+1))
yuw[:-1,1:] = gdepw[:,:,midX] ; zbot = gdepw[-1,:,midX] + 2.*(gdept[-1,:,midX] - gdepw[-1,:,midX])
yuw[:-1,0 ] = gdepw[:,0,midX] ; yuw[-1,1:] = zbot ;  yuw[-1,0] = zbot[0]

xf = np.copy(yuw)*0. ; dx=2*(gphiv[0,0]-gphit[0,0])/1E3
for k in range(nK+1):
    xf[k,1:]=gphiv[:,midX]/1E3
xf[:,0] = xf[:,1] - dx
########################################################
# masking density points
temp = np.zeros((nT,nK,nY))
for t in range(nT):
    temp[t] = toce[t,:,:,midX]
# https://stackoverflow.com/questions/32171917/how-to-copy-a-2d-array-into-a-3rd-dimension-n-times
# indexing with np.newaxis inserts a new 3rd dimension, which we then repeat the
# array along, (you can achieve the same effect by indexing with None, see below)
bigmask = 1-np.repeat(tmask[np.newaxis,:,:,midX],nT,axis=0)
temp = np.ma.array(temp,mask=bigmask)
########################################################
uoce = np.zeros((nT,nK,nY))
for t in range(nT):
    uoce[t] = uu[t,:,:,midX]
# https://stackoverflow.com/questions/32171917/how-to-copy-a-2d-array-into-a-3rd-dimension-n-times
# indexing with np.newaxis inserts a new 3rd dimension, which we then repeat the
# array along, (you can achieve the same effect by indexing with None, see below)
bigmask = 1-np.repeat(umask[np.newaxis,:,:,midX],nT,axis=0)
uoce = np.ma.array(uoce,mask=bigmask)
########################################################
# to have a look at the meshmask
# first bottom cell (masked - task[mbathy,:,:]-> 0)
zht = mbathy*0.
for ii in range(nX):
    for jj in range(nY):
        kk = mbathy[jj,ii]
        zht[jj,ii] = gdepw[kk,jj,ii]
########################################################

tmin = 45.6 ; tmax = 46.4
levelt=np.arange(tmin,tmax,0.05)

vmin = -0.5 ; vmax = 0.5 ; dv = 0.1 # m/s
levelv=np.arange(vmin,vmax,dv)
Ncolor = len(levelv)

ly = np.max((40,1.5*(gphit[midY+1,midX]-gphit[midY,midX])/1E3))
titlezer  = "Cross section"
palette = plt.get_cmap('RdBu_r',Ncolor)
optpcolor = {"vmin":vmin, "vmax":vmax, "cmap" : palette , 'alpha': 0.5}
optcontour = {'levels' : levelv ,"vmin":vmin, "vmax":vmax,'colors':'black', 'linestyles' : "solid"}

opthatch = {'facecolor':'grey', 'alpha' : 0.5, 'interpolate':True, 'step' : 'mid'}
optcontour_t = {'levels' : levelt ,"vmin":tmin, "vmax":tmax,'colors':'black', 'linestyles' : "solid"}

fig, ax = plt.subplots(dpi=200)

c  = ax.contour(xt,yt,temp[0], **optcontour_t)
cl = ax.clabel(c, c.levels, inline=True, fmt = "%.2f", fontsize=8)
cf  = ax.pcolormesh(xf, yuw, uoce[0],**optpcolor)
# cfc = ax.contour(xt,yt,uoce[0], **optcontour)

cbar = plt.colorbar(cf)
cbar.set_label(r"zonal speed (m/s)")

# ... we would expect UW points to be the bottom of the basin
ax.fill_between(gphit[:,0]/1E3, 6000., zht[:,midX], **opthatch)

if timeframe=="12h" :
    titlezer = '%02dd/%02dd \n'%((0.)*tskip/2.,nT*tskip/2.) # because 12h sortie
elif timeframe=="2h":
    titlezer = '%02dh/%02dh (%02dd/%02dd) \n'%((0.)*tskip*2.,nT*tskip*2.,(0.)*tskip/12.,nT*tskip/12.) # because 2h sortie
titlezer += "\nmin = %1.1f m/s   max = %1.1f m/s" % (np.nanmin(uoce[0]), np.nanmax(uoce[0]))
ax.set_title(titlezer, fontsize = 12, y = 1.02)

ax.set_ylim(4800,0)
ax.set_yticks([0,1000,2000,3000,4000])
# ax.set_yticks([0,1000,2000,3000,4000,5000])
ax.set_ylabel("Z (m)")

ax.set_xlim(ly,-ly)
ax.set_xlabel("Y (km)")

ax.xaxis.set_major_locator(MultipleLocator(ly))
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(MultipleLocator(250))
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
plt.tight_layout()

def animate(i):
    """Set the data for the ith iteration of the animation."""
    global c,cl,cf, sfm
    #
    ax.collections = []
    for label in cl:
        label.remove()
    #
    c  = ax.contour(xt,yt,temp[i], **optcontour_t)
    cl = ax.clabel(c, c.levels, inline=True, fmt = "%.2f", fontsize=8)
    cf  = ax.pcolormesh(xf, yuw, uoce[i],**optpcolor)
    # cfc = ax.contour(xt,yt,uoce[i], **optcontour)

    # cf = ax.pcolormesh(xt, yt, sfm[i], **optpcolor)
    #
    ax.fill_between(gphit[:,0]/1E3,6000., zht[:,midX], **opthatch)
    #
    if timeframe=="12h" :
        ptitle = '%02dd/%02dd \n'%((i)*tskip/2.,nT*tskip/2.) # because 12h sortie
    elif timeframe=="2h":
        ptitle = '%02dh/%02dh (%02dd/%02dd) \n'%((i)*tskip*2.,nT*tskip*2.,(i)*tskip/12.,nT*tskip/12.) # because 2h sortie
    ptitle += "\nmin = %1.1f m/s   max = %1.1f m/s" % (np.nanmin(uoce[i]), np.nanmax(uoce[i]))
    sys.stdout.write(u"\u001b[1000D" + "processing movie [%3d/%3d]" % (i+1,nT))
    sys.stdout.flush()
    ax.set_title(ptitle, fontsize = 12, y = 1.02)
    #
    return c,cl,cf

if save:
    anim = animation.FuncAnimation(fig, animate, frames=nT, blit=False, repeat=False)
    writer = animation.writers['ffmpeg'](fps=4)
    anim.save('%s.mp4' % (psave), writer=writer, dpi=200)
    plt.close("all")
    print("\nsaving : %s" % psave)
else:
    anim = animation.FuncAnimation(fig, animate, frames=nT)
    plt.show()
