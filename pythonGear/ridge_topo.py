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
# ok donc r maximale est le delta des h sur la mi hauteur.
# si je résous ma pente

pmm = "mesh_mask.nc"
# pmm = "mesh_mask_dx10_dy5.nc"
#pmm = "mesh_mask_dx10_dy10.nc"
# pmm = "mesh_mask_dx10_dy20.nc"
#
dirr = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref/"
dirr = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref/"
pmm = dirr+"/mesh_mask9.nc"
#
# pmm = "mesh_mask_dx5_dy5.nc"
# rn_dy = 10. ; rn_dx = 10.
#
Lx = 1000E3 ; Ly = 400E3 ; Hz = 2400.
rn_hdepth   = 5500.
rn_h_strait   = 4500. ; rn_h_ridge    = 2800.
rn_std_strait = 30.   ; rn_std_ridge  = 200.   # m  km
rn_width      = 30    ; rn_delta      = 5    # km km
#
zphi0=0.
#
psave = "gtopo.png" ; save = 1
# psave = "gtopo_35.png" ; save = 1

""" Mask
"""
mm = nc4.Dataset(pmm)
tmask = mm.variables['tmask'][0,:,:,:]

nK,nJ,nI = np.shape(tmask)
midJ = nJ//2 # careful, mid on T point are staggered
midI = nI//2 # careful, mid on T point are staggered

tmask = tmask[:,:,midI]
glamt = mm.variables['glamt'][0,:,:] ; glamu = mm.variables['glamu'][0,:,:]
gphit = mm.variables['gphit'][0,:,:] ; gphiu = mm.variables['gphiu'][0,:,:] ; gphiv = mm.variables['gphiv'][0,:,:]
e1t   = mm.variables['e1t']  [0,:,:] ; e2t   = mm.variables['e2t']  [0,:,:]
mbathy= mm.variables['mbathy'][0,:,:] ; gdepw = mm.variables['gdepw_0'][0,:,:,:]

""" Geometrical setting
"""
# to test new geometry
# equals 0 outside and -1 inside
# z1d = 1. + 0.5 * ( np.tanh( (gphit - 1E3*rn_width/2.) / (rn_delta*1E3) )  \
#                  - np.tanh( (gphit + 1E3*rn_width/2.) / (rn_delta*1E3) )  )
# # strait and ridge shapes
# z1h =           z1d * rn_h_ridge   + (1.-z1d) * rn_h_strait
# z1s = 1000. * ( z1d * rn_std_ridge + (1.-z1d) * rn_std_strait ) # in meters as glam
# zht =  rn_hdepth - (rn_hdepth - z1h) * np.exp( - glamt**2 / (2.*z1s*z1s) ) # depuis le fond
# #zht = z1h * np.exp( - glamt**2 / (2.*z1s*z1s) ) # depuis la surface

# to have a look at the meshmask
zht = mbathy*0.
for ii in range(nI):
    for jj in range(nJ):
        zht[jj,ii] =gdepw[mbathy[jj,ii],jj,ii]

zhu = glamt*0. ; zhv = glamt*0.
for i in range(nI-1):
    for j in range(nJ-1):
        zhu[j,i] = 0.5 * (zht[j,i] + zht[j  ,i+1])
        zhv[j,i] = 0.5 * (zht[j,i] + zht[j+1,i  ])

# gradu h / hu ou gradv h / hv glam[y,x]
ru = glamt*0. ; rv = glamt*0. ; rt = glamt*0. ; rmax=glamt*0.
for i in range(nI-1):
    ru[:,i] = np.abs(zht[:,i] - zht[:,i-1]) / (zht[:,i] + zht[:,i-1])

for j in range(nJ-1):
    rv[j,:] = np.abs(zht[j,:] - zht[j-1,:]) / (zht[j,:] + zht[j-1,:])

for i in range(nI-1):
    for j in range(nJ-1):
        rt[j,i] = 0.25 * ( ru[j,i] + ru[j  ,i+1]   \
                         + rv[j,i] + rv[j+1,i  ]   )

for i in range(nI-1):
    for j in range(nJ-1):
        rmax[j,i] = np.max([ru[j,i],ru[j  ,i+1], \
                            rv[j,i],rv[j+1,i  ]] )

gridx = np.zeros((nJ+1,nI+1)) ; gridy = np.zeros((nJ+1,nI+1))
for i in range(0,nI):
    gridx[:-1,i+1] = glamu[:  ,i]
for j in range(1,nJ):
    gridy[j+1,:-1] = gphiv[j,:  ]

gridx[:-1,0] = glamu[:,0]  - e1t[:,0]  ; gridy[0,:-1] = gphiv[0,:] - e2t[0,:]
gridx[ -1,0] = glamu[-1,0] - e1t[-1,0] ; gridy[0, -1] = gphiv[0,-1] - e2t[0,-1]

gridx[-1,1:] = glamu[-1,:]   ; gridy[1:,-1] = gphiv[:,-1]

######

# data = zht
# data = ru
data = rmax
# data = np.max(ru,rv)

######

""" Plot
"""
# Lx = 400E3 ; Ly = 40E3

vmin=rn_h_ridge ;  vmax=rn_hdepth
cticks = np.arange(vmin, vmax, 250)
# palette = plt.get_cmap('YlOrRd',20)
palette = plt.get_cmap('Blues',5)

# https://matplotlib.org/stable/gallery/lines_bars_and_markers/scatter_hist.html#sphx-glr-gallery-lines-bars-and-markers-scatter-hist-py
left, width = 0.1, 0.6
bottom, height = 0.1, 0.6
spacing = 0.02
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]

# start with a square Figure
fig = plt.figure(figsize=(8, 8),dpi=200)

ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

""" SECTIONS X (top panel) """

ax_histx.plot(glamt[0,:], zht[0,:], 'k')
ax_histx.fill_between(glamt[0,:], rn_hdepth, zht[0,:], hatch = "|",
                      facecolor='limegreen', alpha = 0.8,interpolate=True)
ax_histx.plot(glamt[midJ,:], zht[midJ,:], 'k.-')
ax_histx.fill_between(glamt[midJ,:], rn_hdepth, zht[midJ,:], hatch = "|",
                      facecolor='grey', alpha = 1.,interpolate=True)
ax_histx.tick_params(axis="x", labelbottom=False)
ax_histx.set_ylabel("z (m)")
# ax_histx.invert_yaxis()
ax_histx.set_ylim(rn_hdepth,Hz)
ax_histx.axhline(xmin=glamt[midJ,0],xmax=glamt[midJ,-1], y=rn_h_ridge, linestyle = '--', color = 'red', lw = 1.)
ax_histx.axhline(xmin=glamt[midJ,0],xmax=glamt[midJ,-1], y=rn_h_strait, linestyle = '-', color = 'red', lw = 1.)

""" SECTIONS Y (right panel) """
# attention - comme axes renversé, trace mal les marches
# ax_histy.plot(zht[:,midI], gphit[:,midI], color='black', marker='o', markersize = 3)

# ax_histy.plot(zht[:,midI], gphit[:,midI], color='black',drawstyle='steps-mid', marker='o', markersize = 3)
ax_histy.fill_betweenx(gphit[:,midI], rn    hdepth, rn_h_strait,hatch = "-",
                      facecolor='grey', alpha = 0.9,interpolate=True)
ax_histy.fill_betweenx(gphit[:,midI], rn_h_strait, zht[:,midI],hatch = "-",
                      facecolor='limegreen', alpha = 0.9,interpolate=False, step = 'mid')
ax_histy.scatter(zht[:,midI], gphit[:,midI], color='black', marker='o', s = 10)

ax_histy.tick_params(axis="y", labelleft=False)
ax_histy.set_xlabel("z (m)")
# ax_histy.invert_xaxis()
ax_histy.set_xlim(rn_hdepth,Hz)
ax_histy.axvline(ymin=gphit[0,midI],ymax=gphit[-1,midI], x=rn_h_ridge, linestyle = '--', color = 'red', lw = 1.)
ax_histy.axvline(ymin=gphit[0,midI],ymax=gphit[-1,midI], x=rn_h_strait, linestyle = '-', color = 'red', lw = 1.)
# ax_histy.axhline(y=0,linestyle = '--', color = 'limegreen', lw =1.)
ax_histy.axhline(y=0.             , linestyle = '--', color = 'royalblue', lw = 2.)
ax_histy.axhline(y= rn_width*1E3/2, linestyle = '-', color = 'royalblue', lw = 1.)
ax_histy.axhline(y=-rn_width*1E3/2, linestyle = '-', color = 'royalblue', lw = 1.)
ax_histy.axhline(y= gphiv[midJ  ,midI], linestyle = '--', color = 'red', lw = 1.)
ax_histy.axhline(y= gphiv[midJ-1,midI], linestyle = '--', color = 'red', lw = 1.)

""" MAIN FIGURE """
# fig, ax = plt.subplots(dpi=200)
im = ax.pcolormesh(gridx, gridy, data,
                   vmin = 0.,
                   alpha = 0.95,
                   cmap = palette)
c = ax.contour(glamt,gphit,zht, levels = cticks,colors='k')
cbarticks = [0.,np.nanmax(data)]
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
cbaxes = inset_axes(ax, width="30%", height="3%", loc=8)
cbar = plt.colorbar(im,cax=cbaxes, ticks=cbarticks, orientation='horizontal')
cbar.set_ticklabels(["{:.2f}".format(x) for x in cbarticks])
cbar.set_label(r"r max", labelpad = -30)
cbar.ax.tick_params(pad=-30)

c1 = ax.contour(glamt, gphit, zht,
                levels = cticks,
                vmin=vmin,vmax=vmax,
                linewidths =0.4, colors=('k',),linestyles = "solid")
ax.clabel(c1, fmt='%3.0f', colors='k', fontsize=6, inline=True)
ax.set_xlim(-Lx,Lx)
tickszer = [np.float64(i-4)*Lx/4. for i in range(9)]
ax.set_xticks(tickszer)
ax.set_xticklabels(["%d" % (x/1E3) for x in tickszer])
ax.set_xlabel("X (km)")

ax.set_ylim(-Ly,Ly)
tickszer = [np.float64(i-4)*Ly/4. for i in range(9)]
ax.set_yticks(tickszer)
ax.set_yticklabels(["%d" % (y/1E3) for y in tickszer])
ax.set_ylabel("Y (km)")

""" PRETTY """
ax.patch.set_color('0.')
ax.xaxis.set_minor_locator(MultipleLocator(Lx/(4*4)))
ax.yaxis.set_minor_locator(MultipleLocator(Ly/(4*4)))
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5, left = True, right = True)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10, bottom = True, top = True)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
# ax.grid(which='both')

# ax_histx.xaxis.set_minor_locator(AutoMinorLocator())
ax_histx.yaxis.set_minor_locator(MultipleLocator(250.))
ax_histx.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5, left = True, right = True)
ax_histx.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10, bottom = True, top = True)
ax_histx.tick_params(which='minor',length = 4)
ax_histx.tick_params(which='major',length = 6)

ax_histy.xaxis.set_minor_locator(MultipleLocator(250.))
# ax_histy.yaxis.set_minor_locator(AutoMinorLocator())
ax_histy.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5, left = True, right = True)
ax_histy.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10, bottom = True, top = True)
ax_histy.tick_params(which='minor',length = 4)
ax_histy.tick_params(which='major',length = 6)

# plt.tight_layout

if save :
    print("\nsaving : %s" % psave)
    fig.savefig(psave, dpi = 200)
    plt.close()
plt.show()
