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
from matplotlib import colors
import matplotlib.animation as animation
# from pylab import *
# import cmocean


""" *****************************************************************
"""

dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref"
dirm = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref"

# pdt = "/RIDGE_sco_1_12h_grid_T.nc"
# pdu = "/RIDGE_sco_1_12h_grid_U.nc"
# pdv = "/RIDGE_sco_1_12h_grid_V.nc"
pdw = "/RIDGE_sco_1_12h_grid_W.nc"
pmm = "/mesh_mask.nc"

save = 1 ; psave = "ridgekz" ; film = 1

########################################################
vmin = 1.e-5 ; vmax = 1.e-3 ; tskip=4
Ncolor=8
dv=(vmax-vmin)/Ncolor
"""
    Plot du test overflow - zps : z partial cells
    34 : time step
    101 : depth (z)
    3 : y dimension (1 cell embraced by two boundary layers)
    202 : x dimension
"""
dtw = nc4.Dataset(dir+pdw)
# dt  = nc4.Dataset(dir+pdt)
mm = nc4.Dataset(dirm+pmm)

kz   = dtw.variables['avt'][::tskip,:,:,:]   # volumic transport
# toce = dt.variables['toce'] [::tskip,:,:,:]
nT,nK,nY,nX = np.shape(kz)
# nT=1

tmask = mm.variables['tmask'][0][:,:,:]
glamt = mm.variables['glamt'][0] ; gphit = mm.variables['gphit'][0]
glamu = mm.variables['glamu'][0]
umask = mm.variables['umask'][0] ; vmask = mm.variables['vmask'][0]
gdept = mm.variables['gdept_0'][0] ; gdepw = mm.variables['gdepw_0'][0]
e3t = mm.variables['e3t_0'][0]

print("domain size is (x,y) %dx%d with %d k-levels" % (nX,nY,nK))
midY = nY//2# np.where(np.abs(gphit[:,0])<=1E3)[0][0]
# attention en BVP le milieu est un V point
midX = np.where(np.abs(glamt[0,:])<=1E3)[0][0]

#######################################################
# kz are W points
yw = gdepw[:,midY,:-1]
for i in range(nX-1):
    yw[:,i]=0.5*( gdepw[:,midY,i] + gdepw[:,midY,i+1] )
xt = yw*0.
for k in range(nK):
    xt[k,:]=glamt[midY,:-1]/1E3

# ... surrounded by U points
dz = e3t[0,midY,:]
yt = np.zeros((nK+1,nX))
yt[0,:] = 0. ; yt[1:,:] = gdept[:,midY,:]
xu = yt*0. ; x1t = glamu[  midY,:]/1E3
for k in range(nK+1):
    xu[k,:]=x1t

########################################################
avt = np.zeros((nT,nK,nX))
for t in range(nT):
    avt[t] = kz[t,:,midY,:]
# https://stackoverflow.com/questions/32171917/how-to-copy-a-2d-array-into-a-3rd-dimension-n-times
# indexing with np.newaxis inserts a new 3rd dimension, which we then repeat the
# array along, (you can achieve the same effect by indexing with None, see below)
bigmask = 1-np.repeat(tmask[np.newaxis,:,midY,:],nT,axis=0)
avt = np.ma.array(avt,mask=bigmask)
avt[avt==0.] = 1.e-5
########################################################

titlezer  = "kz"
palette = plt.get_cmap('RdBu_r',Ncolor)
optpcolor = {"vmin":vmin, "vmax":vmax, "cmap" : palette}
opthatch = {'facecolor':'grey', 'alpha' : 0.5, 'interpolate':True}
def ridge(x):
    return(5500-(5500-3000)*np.exp(-x**2/(2.*200**2)))
# def strait(x):
#     return(5500-(5500-4500)*np.exp(-x**2/(2.*30**2)))
def strait(x):
    return(5500-(5500-4500)*np.exp(-x**2/(2.*150**2)))

fig, ax = plt.subplots(dpi=200)
# cf = ax.pcolormesh(xu, yt, avt[0],**optpcolor)
# c  = ax.contour(xt,yw,avt[0], levels = Ncolor//2, colors='k')
from matplotlib.colors import LogNorm
Z = avt[-1]
cf = ax.pcolormesh(xu,yt,Z, **optpcolor, norm=colors.LogNorm())
cbar = plt.colorbar(cf)
cbar.set_label(r"log(kz) (m2/s)")
# ax.fill_between(xu[0,:], strait(xw[0,:]), ridge(xw[0,:]),**opthatch)
# ... we would expect UW points to be the bottom of the basin
# ax.fill_between(xt[0,:], strait(xt[0,:]), ridge(xt[0,:]),**opthatch)

titlezer = '%02d/%02d \n'%((0.)*tskip/2.,nT*tskip/2.) # because 12h sortie
titlezer += "\nmin = %1.1e m2/s   max = %1.2e m2/s" % (np.nanmin(Z), np.nanmax(Z))
ax.set_title(titlezer, fontsize = 12, y = 1.02)

ax.set_ylim(5500,0)
ax.set_yticks([0,1000,2000,3000,4000,5000])
# ax.set_yticks([0,1000,2000,3000,4000,5000])
ax.set_ylabel("Z (m)")

ax.set_xlim(-1000,1000)
ax.set_xlabel("X (km)")

ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.yaxis.set_minor_locator(MultipleLocator(250))
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
plt.tight_layout()

# def animate(i):
#     """Set the data for the ith iteration of the animation."""
#     global c,cf, sfm
#
#     ax.collections = []
#     cf = ax.pcolormesh(xu, yt, avt[i], **optpcolor)
#     c = ax.contour(xt,yw,avt[i], levels = Ncolor//2,colors='k')
#     #
#     ax.fill_between(xt[0,:], strait(xt[0,:]), ridge(xt[0,:]),**opthatch)
#     #
#     ptitle = '%02dd/%02dd \n'%((i+1)*tskip/2.,nT*tskip/2.) # because 12h sortie
#     ptitle += "min = %2.2f Sv   max = %2.2f Sv   " % (np.nanmin(sfm[i]), np.nanmax(sfm[i]))
#     sys.stdout.write(u"\u001b[1000D" + "processing movie [%3d/%3d]" % (i+1,nT))
#     sys.stdout.flush()
#     ax.set_title(ptitle, fontsize = 12, y = 1.02)
#     #
#     return c, cf
#
# if save:
#     anim = animation.FuncAnimation(fig, animate, frames=nT, blit=False, repeat=False)
#     writer = animation.writers['ffmpeg'](fps=4)
#     anim.save('%s.mp4' % (psave), writer=writer, dpi=200)
#     plt.close("all")
#     print("\nsaving : %s" % psave)
# else:
#     anim = animation.FuncAnimation(fig, animate, frames=nT)
#     plt.show()
if save :
    print("saving : %s" % psave)
    fig.savefig(psave)
    plt.close()
else :
    plt.show()
