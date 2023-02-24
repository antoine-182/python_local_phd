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

dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref/"
dirm = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref/1ubs"

# dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_bvp"
# dirm = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_bvp"


timeframe = "2h"
pdt = "/RIDGE_sco_1_%s_grid_T.nc" % (timeframe)
pdu = "/RIDGE_sco_1_%s_grid_U.nc" % (timeframe)
pdv = "/RIDGE_sco_1_%s_grid_V.nc" % (timeframe)
pdw = "/RIDGE_sco_1_%s_grid_W.nc" % (timeframe)
# pdw = "/RIDGE_ref_1_12h_grid_W.nc"
pmm = "/mesh_mask.nc"

save = 1 ; psave = "ridge" ; film = 1 ; Nctd = "F"
tskip=100

########################################################
"""
    Plot du test overflow - zps : z partial cells
    34 : time step
    101 : depth (z)
    3 : y dimension (1 cell embraced by two boundary layers)
    202 : x dimension
"""
dtu = nc4.Dataset(dir+pdu)
# dtt = nc4.Dataset(dir+pdt)
dtw = nc4.Dataset(dir+pdw)
mm = nc4.Dataset(dirm+pmm)

# toce = dtt.variables['toce'][:ntend:tskip,:,:,:]   # volumic transport
# kz   = dtw.variables['avt'][:ntend:tskip,:,:,:]   # vertical mixing
uu   = dtu.variables['u_vol'][::tskip,:,:,:]   # volumic transport

nT,nK,nY,nX = np.shape(uu)
# nT=1

tmask = mm.variables['tmask'][0][:,:,:]
# barotropic streamfunction
glamt = mm.variables['glamt'][0] ; gphit = mm.variables['gphit'][0]
# density
glamf = mm.variables['glamf'][0] ; gphif = mm.variables['gphif'][0]
glamu = mm.variables['glamu'][0] ; gphiv = mm.variables['gphiv'][0]
umask = mm.variables['umask'][0] ; gdepw = mm.variables['gdepw_0'][0]
mbathy= mm.variables['mbathy'][0,:,:]

print("domain size is (x,y) %dx%d with %d k-levels" % (nX,nY,nK))
if Nctd == "F" :
    midY = np.where(np.abs(gphiv[:,0])<=1E3)[0][0] # same T point
elif Nctd == "U" :
    midY = np.where(np.abs(gphit[:,0])<=1E3)[0][0] # same T point

midX = np.where(np.abs(glamt[0,:])<=1E3)[0][0]

########################################################
# to have a look at the meshmask
# first bottom cell (masked - task[mbathy,:,:]-> 0)
zht = mbathy*0.
for ii in range(nX):
    for jj in range(nY):
        kk = mbathy[jj,ii]
        zht[jj,ii] = gdepw[kk,jj,ii]
zht = np.ma.array(zht,mask=1-tmask[0])

#######################################################
# barotropic streamfunction nodes are FW points # xu and yv
yv = np.copy(gphiv[:,:])/1E3
xu = np.copy(glamu[:,:])/1E3

# ... surrounded by T points
yt = np.zeros((nY+1,nX+1))
for ii in range(nX):
    yt[:-1,ii]=gphit[:,ii]
#uppest
dy = 2.*(gphiv[-1,:]-gphit[-1,:])
yt[-1,:-1] = gphit[-1,:]+dy
#eastest
yt[:,-1] = yt[:,-2]

xt = np.zeros((nY+1,nX+1))
for jj in range(nY):
    xt[jj,:-1]=glamt[jj,:]
#uppest
dx = 2.*(glamu[:,-1]-glamt[:,-1])
xt[:-1,-1] = glamt[:,-1]+dx
#eastest
xt[-1,:] = xt[-2,:]
#norm
yt/=1E3 ; xt /=1E3
########################################################
sfb = np.zeros((nT,nY,nX)) # u(:,x,y:)
""" streamfunction """
bigmask = 1-np.repeat(umask[np.newaxis,:,:,:],nT,axis=0)
uu = np.ma.array(uu,mask=bigmask)
for jj in range(nY):
    for ii in range(nX):
        sys.stdout.write(u"\u001b[1000D" + "processing psi [%3d/%3d]" % (1+ii+jj*nX,nY*nX))
        sys.stdout.flush()
        sfb[:,jj,ii] += - np.sum(uu[:,:,:jj+1,ii],axis=(1,2) )/1E6

# for t in range(nT):
#     sys.stdout.write(u"\u001b[1000D" + "processing psi [%3d/%3d]" % (t+1,nT))
#     sys.stdout.flush()
#     uum  = np.ma.masked_where(umask==0,uu[t])
#     for jj in range(nY):
#         for ii in range(nX):
#             # k=0 surface - k=nK-1 last dot (masked)
#             # there are as many u than psi
#             sfb[t,jj,ii] += - np.sum(uum[:,:jj+1,ii])/1E6

bigmask = 1-np.repeat(umask[np.newaxis,0,:,:],nT,axis=0)
sfb = np.ma.array(sfb,mask=bigmask)
########################################################
# density points are T points # xt and yt
# yt = np.copy(gphit[:,:])/1E3
# xt = np.copy(glamt[:,:])/1E3
#
# # ... surrounded by F points
# yf = np.zeros((nY+1,nX+1))
# for ii in range(nX):
#     yf[1:,ii]=gphif[:,ii]
# #lowest
# dy = 2.*(gphiv[0,:]-gphit[0,:])
# yf[0,1:] = gphif[0,:]-dy
# #western
# yf[:,0] = yf[:,1]
#
# xf = np.zeros((nY+1,nX+1))
# for jj in range(nY):
#     xf[jj,1:]=glamf[jj,:]
# #western
# dx = 2.*(glamf[:,0]-glamt[:,0])
# xf[1:,0] = glamt[:,0]-dx
# #eastest
# xf[0,:] = xf[1,:]
# #norm
# yf/=1E3 ; xf/=1E3
########################################################
# # bottom temperature
# temp = np.zeros((nT,nY,nX)) # u(:,x,y:)
# for ii in range(nX):
#     for jj in range(nY):
#         temp[:,jj,ii] = toce[:,mbathy[jj,ii]-1,jj,ii]
# temp = np.ma.array(temp,mask=bigmask)

# anomalies
# temp = np.zeros((nT,nY,nX))
# tleft = toce[0,:,midY,1] ; tright = toce[0,:,midY,-2]
# for ii in range(nX):
#     for jj in range(nY):
#     #
#     if ii<=midX :    # left anomaly
#         for jj in range(nY):
#             temp[:,jj,ii] = np.max(tleft*tmask[:,jj,ii] - toce[:,:,jj,ii])
#     elif ii>midX :    # right anomaly
#         for jj in range(nY):
#             temp[:,jj,ii] = np.max(tright*tmask[:,jj,ii] - toce[:,:,jj,ii])
# tinit = toce[0,:,:,:] * tmask
# temp[:,:,:] = np.max(tinit - toce[:,:,:,:], axis=(1,))
# temp = np.ma.array(temp,mask=bigmask)

# data = np.zeros((nT,nY,nX))
# data = np.max(kz[:,:,:,:], axis=(1,))
# data = np.ma.array(data,mask=bigmask)

########################################################
cmin = -6. ; cmax = 6. ; Nc=12             # contour stream
levelc = np.linspace(cmin,cmax, Nc//2+1)

pmin = -6. ; pmax = 6. ; Np=9             # pcolor stream
levelp = np.linspace(pmin,pmax, Np//2+1)

# pmin = 45.7; pmax = 46.4 ; Np = 14        # pcolor bottom temp
# levelp = np.linspace(pmin,pmax, Np//2+1)

# pmin = 1e-5; pmax = 100 ; Np = 14        # pcolor anomaly
# levelp = np.logspace(pmin,pmax, Np)

# streamfunction
datap  = sfb ; datac = sfb
xcorner = xt ; xdata = xu
ycorner = yt ; ydata = yv
titlebar = r"Streamfunction $\psi$ (Sv)"
palette = plt.get_cmap('RdBu_r',Np)

# temp + streamfunction
# datap  = temp ; datac = sfb
# xcorner = xf  ; xdata = xt
# ycorner = yf  ; ydata = yt
# palette = plt.get_cmap('RdBu',Np)
# titlebar = r"Bottom density anomaly (kg/m3)"

# kz + streamfunction
# datap  = data ; datac = sfb
# xcorner = xf  ; xdata = xt
# ycorner = yf  ; ydata = yt
# palette = plt.get_cmap('YlGnBu',Np)
# titlebar = r"Vertical diffusivity"
# import matplotlib.colors as colors
########################################################

titlezer  = "Barotropic Streamfunction"
optpcolor = {"vmin":pmin, "vmax":pmax, "cmap" : palette, "alpha" : 0.5}
# optpcolor = {"cmap" : palette}

opthatch   = {'facecolor':'grey', 'alpha' : 0.5, 'interpolate':True, 'step' : 'mid'}
optclabel  = {'inline' : 'True', 'fmt'  :  "%1.1f", 'fontsize' : 8}
optzlabel  = {'inline' : 'True', 'fmt'  :  "%dm", 'fontsize' : 6}
optcontour = {'levels' : levelc, 'colors':'black', 'linestyles' : "solid", }
# optcontour = {'colors':'black', 'linestyles' : "solid"}
cticks = np.arange(2500, 5500, 500)

fig, ax = plt.subplots(dpi=200)
zc = ax.contour(glamt/1E3,gphit/1E3,zht, levels = cticks,colors='grey', linewidths = 0.5)
# ax.clabel(zc, zc.levels, **optzlabel)

c  = ax.contour( xdata, ydata, datac[0], **optcontour)
cl = ax.clabel(c, c.levels, **optclabel)

# cf = ax.pcolormesh(xcorner,ycorner,datap[0],**optpcolor, norm=colors.LogNorm(vmin=datap.min(), vmax=datap.max()))
cf = ax.pcolormesh(xcorner,ycorner,datap[0],**optpcolor)

cbar = plt.colorbar(cf)
cbar.set_label("%s" % (titlebar))
# ax.fill_between(xu[0,:], strait(xw[0,:]), ridge(xw[0,:]),**opthatch)
# ... we would expect UW points to be the bottom of the basin
# ax.fill_between(glamt[midY,:]/1E3, gphit[:,midX]/1E3, zht[midY,:], **opthatch)
if timeframe=="12h" :
    titlezer = '%02dd/%02dd \n'%((0.)*tskip/2.,nT*tskip/2.) # because 12h sortie
elif timeframe=="2h":
    titlezer = '%02dh/%02dh (%02dd/%02dd) \n'%((0.)*tskip*2.,nT*tskip*2.,(0.)*tskip/12.,nT*tskip/12.) # because 2h sortie
titlezer += "\nmin = %1.1f Sv   max = %1.1f Sv" % (np.nanmin(datac[0]), np.nanmax(datac[0]))
ax.set_title(titlezer, fontsize = 12, y = 1.02)

ax.set_ylim(-400,400)
# ax.set_yticks([2000,3000,4000,5000])
# ax.set_yticks([0,1000,2000,3000,4000,5000])
ax.set_ylabel("Y (km)")

ax.set_xlim(-1000,1000)
ax.set_xlabel("X (km)")

# ax.xaxis.set_major_locator(MultipleLocator(500))
# ax.xaxis.set_minor_locator(MultipleLocator(100))
# ax.yaxis.set_minor_locator(MultipleLocator(250))
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
plt.tight_layout()

def animate(i):
    """Set the data for the ith iteration of the animation."""
    global c,cf,cl

    ax.collections = []
    for label in cl:
        label.remove()

    ax.contour(glamt/1E3,gphit/1E3,zht, levels = cticks,colors='grey', linewidths = 0.5)
    c  = ax.contour( xdata, ydata, datac[i], **optcontour)
    cl = ax.clabel(c, c.levels, **optclabel)
    cf = ax.pcolormesh(xcorner,ycorner,datap[i],**optpcolor )
    #
    if timeframe=="12h" :
        ptitle = '%02dd/%02dd \n'%((i)*tskip/2.,nT*tskip/2.) # because 12h sortie
    elif timeframe=="2h":
        ptitle = '%02dh/%02dh (%02dd/%02dd) \n'%((i)*tskip*2.,nT*tskip*2.,(i)*tskip/12.,nT*tskip/12.) # because 2h sortie
    ptitle += "min = %2.2f Sv   max = %2.2f Sv   " % (np.nanmin(datac[i]), np.nanmax(datac[i]))
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
