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

dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_bvp"
dirm = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_bvp"

# dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_bvp"
# dirm = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_bvp"


pdt = "/RIDGE_sco_1_12h_grid_T.nc"
pmm = "/mesh_mask.nc"

save = 1 ; psave = "ridge" ; film = 1

dt  = nc4.Dataset(dir+pdt)
mm = nc4.Dataset(dir+pmm)

tmask = mm.variables['tmask'][0][:,:,:]
nK,nY,nX = np.shape(tmask)
print("domain size is (x,y) %dx%d with %d k-levels" % (nX,nY,nK))


glamt = mm.variables['glamt'][0] ; gphit = mm.variables['gphit'][0]
glamu = mm.variables['glamu'][0]
umask = mm.variables['umask'][0] ; vmask = mm.variables['vmask'][0]
gdept = mm.variables['gdept_0'][0] ; gdepw = mm.variables['gdepw_0'][0]
e3w = mm.variables['e3w_0'][0]
mbathy= mm.variables['mbathy'][0,:,:]

dy = gphit[1,0] - gphit[0,0]
midY = np.where(np.abs(gphit[:,0])<=dy)[0][0]
midX = np.where(np.abs(glamt[0,:])<=1E3)[0][0]
midX = np.where(np.abs(glamt[0,:])<=1E3)[0][0]

try:
    rpot = dt.variables['rpot'][0,:,:,midX]
except:
    rpot = tmask

#######################################################
# porosity points are T points
yt = np.zeros((nK,nY))
yt[:,:] = gdept[:,:,midX]
xt = yt*0. ; x1t = gphit[:,midX]/1E3
for k in range(nK):
    xt[k,:]=x1t

# ... surrounded by WF points
# yw = np.zeros((nK+1,nX+1))
# yw[:-1,1:] = gdepw[:,midY,:] ; zbot = gdepw[-1,midY,:] + 2.*(gdept[-1,midY,:] - gdepw[-1,midY,:])
# yw[:-1,0 ] = gdepw[:,midY,0] ; yw[-1,1:] = zbot ;  yw[-1,0] = zbot[0]
# yuw= yw*0.
# for i in range(0,nX):
#     yuw[:,i]=0.5*( yw[:,i] + yw[:,i+1] )
#
# xu = yw*0. ; dx=2*(glamu[0,0]-glamt[0,0])/1E3
# for k in range(nK+1):
#     xu[k,1:]=glamu[midY,:]/1E3
# xu[:,0] = xu[:,1] - dx

########################################################
# to have a look at the meshmask
# first bottom cell (masked - task[mbathy,:,:]-> 0)
zht = mbathy*0.
for ii in range(nX):
    for jj in range(nY):
        kk = mbathy[jj,ii]
        zht[jj,ii] = gdepw[kk,jj,ii]

""" Map
    A good manner to test if the grid match the correct data, is to
    plot without masking the datafield.
    vmin and vmax are necessary to turn blank 0 porosity
"""
palette = plt.get_cmap('Blues')
# fig, ax = plt.subplots(figsize = [12, 8])
fig, ax = plt.subplots(dpi=200)

# rpot = np.ma.masked_where(tmask==0,rpot)

im = plt.pcolormesh(xt, yt, rpot,
                   vmin=0,vmax=1,
                   cmap = palette)
ax.plot(gphit/1E3,zht, color='tab:red', marker='o',
        linewidth = 0.8, markersize = 1.5)
# ax.fill_between(gphit[:,0]/1E3, 6000., zht[:,midX], **opthatch)

# cbar = plt.colorbar(im)

ax.patch.set_color('peru')
ax.set_ylabel("Z (100m)")
ax.set_xlabel("Y (km)")

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
plt.tight_layout()

plt.tight_layout()
# ax.set_aspect(aspect=1) # data coordinate 'equal'

##################################################
