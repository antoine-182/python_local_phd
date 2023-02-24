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

dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref/ref3ubs"
dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref/ref3ubs"

# dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_bvp"
# dirm = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_bvp"


timeframe = "2h"
pdt = "/RIDGE_sco_1_%s_grid_T.nc" % (timeframe)
pdu = "/RIDGE_sco_1_%s_grid_U.nc" % (timeframe)
pdv = "/RIDGE_sco_1_%s_grid_V.nc" % (timeframe)
pdw = "/RIDGE_sco_1_%s_grid_W.nc" % (timeframe)
# pdw = "/RIDGE_ref_1_12h_grid_W.nc"
pmm = "/mesh_mask.nc"

save = 1 ; psave = "ridge" ; film = 1

########################################################
vmin = 45.7 ; vmax = 46.3 ; tskip=3
Ncolor=12
dv=(vmax-vmin)/Ncolor
"""
    Plot du test overflow - zps : z partial cells
    34 : time step
    101 : depth (z)
    3 : y dimension (1 cell embraced by two boundary layers)
    202 : x dimension
"""
# dtu = nc4.Dataset(dir+pdu)
dt  = nc4.Dataset(dir+pdt)
mm = nc4.Dataset(dir+pmm)

tmask = mm.variables['tmask'][0][:,:,:]
nK,nY,nX = np.shape(tmask)

glamt = mm.variables['glamt'][0] ; gphit = mm.variables['gphit'][0]
glamu = mm.variables['glamu'][0] ; gphiv = mm.variables['gphiv'][0]
umask = mm.variables['umask'][0] ; vmask = mm.variables['vmask'][0]
gdept = mm.variables['gdept_0'][0] ; gdepw = mm.variables['gdepw_0'][0]
e3w = mm.variables['e3w_0'][0]
mbathy= mm.variables['mbathy'][0,:,:]

# uu   = dtu.variables['u_vol'][::tskip,:,:,:]   # volumic transport
toce = dt.variables['toce'][::tskip,:,:,:]
nT,_,_,_ = np.shape(toce)
# nT=1

print("domain size is (x,y) %dx%d with %d k-levels" % (nX,nY,nK))
midY = np.where(np.abs(gphiv[:,0])<=1E3)[0][0]+1
# attention en BVP le milieu est un V point
midX = np.where(np.abs(glamt[0,:])<=1E3)[0][0]

#######################################################
# s coordinate is painful
# density points are T points
yt = np.zeros((nK,nX))
yt[:,:] = gdept[:,midY,:]
xt = yt*0. ; x1t = glamt[midY,:]/1E3
for k in range(nK):
    xt[k,:]=x1t

# ... surrounded by UW points
yw = np.zeros((nK+1,nX+1))
yw[:-1,1:] = gdepw[:,midY,:] ; zbot = gdepw[-1,midY,:] + 2.*(gdept[-1,midY,:] - gdepw[-1,midY,:])
yw[:-1,0 ] = gdepw[:,midY,0] ; yw[-1,1:] = zbot ;  yw[-1,0] = zbot[0]
yuw= yw*0.
for i in range(0,nX):
    yuw[:,i]=0.5*( yw[:,i] + yw[:,i+1] )

xu = yw*0. ; dx=2*(glamu[0,0]-glamt[0,0])/1E3
for k in range(nK+1):
    xu[k,1:]=glamu[midY,:]/1E3
xu[:,0] = xu[:,1] - dx

########################################################
# masking density points
temp = np.zeros((nT,nK,nX)) ; bigmask = 1-np.repeat(tmask[np.newaxis,:,:,:],nT,axis=0)
toce = np.ma.array(toce,mask=bigmask)
for t in range(nT):
    # temp[t] = toce[t,:,midY,:]
    temp[t] = np.ma.mean(toce[t,:,:,:],axis=1)
# https://stackoverflow.com/questions/32171917/how-to-copy-a-2d-array-into-a-3rd-dimension-n-times
# indexing with np.newaxis inserts a new 3rd dimension, which we then repeat the
# array along, (you can achieve the same effect by indexing with None, see below)
bigmask = 1-np.repeat(tmask[np.newaxis,:,midY,:],nT,axis=0)
temp = np.ma.array(temp,mask=bigmask)
########################################################
# to have a look at the meshmask
# first bottom cell (masked - task[mbathy,:,:]-> 0)
zht = mbathy*0.
for ii in range(nX):
    for jj in range(nY):
        kk = mbathy[jj,ii]
        zht[jj,ii] =gdepw[kk,jj,ii]

########################################################

titlezer  = "Density classes"
palette = plt.get_cmap('RdBu_r',Ncolor)
optpcolor = {"vmin":vmin, "vmax":vmax, "cmap" : palette}
# optpcolor = {"cmap" : palette}
opthatch = {'facecolor':'grey', 'alpha' : 0.5, 'interpolate':True, 'step' : 'mid'}
levelsc = np.linspace(vmin,vmax, Ncolor+1)

fig, ax = plt.subplots(dpi=200)
cf = ax.pcolormesh(xu, yuw, temp[0],**optpcolor)
c  = ax.contour(xt,yt,temp[0], levels = levelsc, colors='k')

cbar = plt.colorbar(cf)
cbar.set_label(r"Density classes (kg/m3)")
# ax.fill_between(xu[0,:], strait(xw[0,:]), ridge(xw[0,:]),**opthatch)
# ... we would expect UW points to be the bottom of the basin
ax.fill_between(xt[0,:], zht[0,:], zht[midY,:], **opthatch)
# ax.fill_between(xt[0,:], 5500, ridge(xt[0,:]),**opthatch)

if timeframe=="12h" :
    titlezer = '%02dd/%02dd \n'%((0.)*tskip/2.,nT*tskip/2.) # because 12h sortie
elif timeframe=="2h":
    titlezer = '%02dh/%02dh (%02dd/%02dd) \n'%((0.)*tskip*2.,nT*tskip*2.,(0.)*tskip/12.,nT*tskip/12.) # because 2h sortie
titlezer += "\nmin = %1.1f   max = %1.1f" % (np.nanmin(temp[0]), np.nanmax(temp[0]))
ax.set_title(titlezer, fontsize = 12, y = 1.02)

ax.set_ylim(5500,1500)
ax.set_yticks([2000,3000,4000,5000])
# ax.set_yticks([0,1000,2000,3000,4000,5000])
ax.set_ylabel("Z (m)")

ax.set_xlim(-1000,1000)
ax.set_xlabel("X (km)")

# ax.patch.set_color('0.')
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.yaxis.set_minor_locator(MultipleLocator(250))
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
plt.tight_layout()

def init():
    return(ax)

def animate(i):
    """Set the data for the ith iteration of the animation."""
    global c,cf, temp

    ax.collections = []
    cf = ax.pcolormesh(xu, yuw, temp[i], **optpcolor)
    c = ax.contour(xt,yt,temp[i], levels = levelsc,colors='k')
    #
    # ax.fill_between(xt[0,:], strait(xt[0,:]), ridge(xt[0,:]),**opthatch)
    ax.fill_between(xt[0,:], zht[0,:], zht[midY,:], **opthatch)
    #
    if timeframe=="12h" :
        ptitle = '%02dd/%02dd \n'%((i)*tskip/2.,(nT-1)*tskip/2.) # because 12h sortie
    elif timeframe=="2h":
        ptitle = '%02dh/%02dh (%02dd/%02dd) \n'%((i)*tskip*2.,(nT-1)*tskip*2.,(i)*tskip/12.,(nT-1)*tskip/12.) # because 2h sortie
    ptitle += "min = %2.2f   max = %2.2f   " % (np.nanmin(temp[i]), np.nanmax(temp[i]))
    sys.stdout.write(u"\u001b[1000D" + "processing movie [%3d/%3d]" % (i,nT-1))
    sys.stdout.flush()
    ax.set_title(ptitle, fontsize = 12, y = 1.02)
    #
    return c, cf

if save:
    anim = animation.FuncAnimation(fig, animate, frames=nT, init_func=init, blit=False, repeat=False)
    writer = animation.writers['ffmpeg'](fps=4)
    anim.save('%s.mp4' % (psave), writer=writer, dpi=200)
    plt.close("all")
    print("\nsaving : %s" % psave)
else:
    anim = animation.FuncAnimation(fig, animate, frames=nT)
    plt.show()
