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
import types

""" *****************************************************************
"""

# ########################################################
#
dir  = "/Users/gm/Documents/nemo/release-4.0.1/tests/my_overflow/EXP_gre"
dirm = "/Users/gm/Documents/nemo/release-4.0.1/tests/my_overflow/EXP_gre"

pdt = dir +"/OVF_grid_T.nc"
pdu = dir +"/OVF_grid_U.nc"
pdw = dir +"/OVF_grid_W.nc"

pmm = dirm+"/mesh_mask.nc"

save = 0 ; psave = "overflow" ; film = 1
nt = -1 ; tskip=1

vmin = -1.5 ; vmax = 1.5 ; timeframe = "30min"
Ncolor=8
dv=(vmax-vmin)/Ncolor


"""
    Plot du test overflow - zps : z partial cells
    34 : time step
    101 : depth (z)
    3 : y dimension (1 cell embraced by two boundary layers)
    202 : x dimension
"""
dtu = nc4.Dataset(pdu)

uu  = dtu.variables['u_vol'][:nt:tskip,:,:,:]   # volumic transport
nT,nK,nY,nX = np.shape(uu)
mm = nc4.Dataset(pmm)

tmask = mm.variables['tmask'][0][:,:,:]
glamt = mm.variables['glamt'][0] ; gphit = mm.variables['gphit'][0]
glamu = mm.variables['glamu'][0]
umask = mm.variables['umask'][0] ; vmask = mm.variables['vmask'][0]
gdept = mm.variables['gdept_0'][0] ; gdepw = mm.variables['gdepw_0'][0]
e3w = mm.variables['e3w_0'][0]
mbathy= mm.variables['mbathy'][0,:,:]


print("domain size is (x,y) %dx%d with %d k-levels" % (nX,nY,nK))
midY = 1
# middle = 3

#######################################################

# streamfunction are UW points
yw = np.copy(gdepw[:,midY,:])
for i in range(nX-1):
    yw[:,i]=0.5*( gdepw[:,midY,i] + gdepw[:,midY,i+1] )
xu = yw*0.
for k in range(nK):
    xu[k,:]=glamu[midY,:]

# ... surrounded by T points
dz = np.copy(e3w[0,midY,:])
yt = np.zeros((nK+1,nX))
yt[0,:] = 0. ; yt[1:,:] = gdept[:,midY,:]
xt = yt*0. ; x1t = glamt[  midY,:]
for k in range(nK+1):
    xt[k,:]=x1t

########################################################
# to have a look at the meshmask
# first bottom cell (masked - task[mbathy,:,:]-> 0)
zht = mbathy*0.
for ii in range(nX):
    for jj in range(nY):
        kk = mbathy[jj,ii]
        zht[jj,ii] = gdepw[kk,jj,ii]

########################################################
sfm = np.zeros((nT,nK,nX))
bigmask = 1-np.repeat(umask[np.newaxis,:,:,:],nT,axis=0)
uu = np.ma.array(uu,mask=bigmask)
for k in range(nK):
    sys.stdout.write(u"\u001b[1000D" + "processing psi [%3d/%3d]" % (k,nK))
    sys.stdout.flush()
    sfm[:,k,:] += - np.sum(uu[:,k:,:,:],axis=(1,2))/1E6

bigmask = 1-np.repeat(umask[np.newaxis,:,1,:],nT,axis=0)
sfm = np.ma.array(sfm,mask=bigmask)

# titlezer  = '%s\n'%(pdu)
# data = sfm[nt,:,:]
########################################################
def pretty(ax,ptitle):
    # ax.plot(xt[slice(10)].T,yt[slice(10)].T*20., 'k-', lw=0.5)
    # https://stackoverflow.com/questions/31877353/overlay-an-image-segmentation-with-numpy-and-matplotlib
    # plt.pcolormesh(gridx, gridk, sill ,alpha=0.9,zorder = 0.)
    # plt.pcolormesh(gridx, gridk, gauss,alpha=0.5,zorder = 0.)
    ax.set_title(ptitle, fontsize = 18, y = 1.02)
    ax.set_ylim(2040,0)
    ax.set_yticks([0,500,1000,1500,2000])
    ax.set_yticklabels(["0","500","1000","1500","2000"])
    ax.set_ylabel("depth (m)")
    ax.set_xlim(0,200)
    ax.set_xticks([0,50, 100, 150, 200])
    ax.set_xticklabels(["0","50","100","150","200"])
    ax.set_xlabel("length (km)")
    ax.patch.set_color('0.8')
    # ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
    ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
    ax.tick_params(which='minor',length = 4)
    ax.tick_params(which='major',length = 6)
    return(ax)

palette = plt.get_cmap('RdBu_r',Ncolor)
optpcolor = {"vmin":vmin, "vmax":vmax, "cmap" : palette, "alpha" : 0.8}
levelsc = np.linspace(vmin,vmax, Ncolor+1)
optcontour = {'levels' : levelsc, 'colors':'black', 'linestyles' : "solid"}

""" figure bordelum """
fig, ax = plt.subplots(figsize = [12, 8])

im = plt.pcolormesh(xt, yt, sfm[0], **optpcolor)
c1 = plt.contour(xu, yw, sfm[0],**optcontour)
titlezer = "min = %2.2f Sv   max = %2.2f Sv   " % (np.min(sfm[0]), np.max(sfm[0]))
cbar = plt.colorbar(im)
pretty(ax,titlezer)
# il faudrait superposer le mask T avec la bonne grille pour vraiment voir ou passe la bathymétry

##################################################

if film :
    """ movie thing """
    def animate(i):
        """Set the data for the ith iteration of the animation."""
        global c,cf, sfm

        ax.collections = []
        c = ax.contour(xu,yw,sfm[i], **optcontour)
        cf = ax.pcolormesh(xt, yt, sfm[i], **optpcolor)
        #
        if timeframe=="30min":
            ptitle = '%02.1fh/%02.1fh \n'%((i)*tskip*0.5,(nT-1)*tskip*0.5) # because 2h sortie
        ptitle += "min = %2.2f Sv   max = %2.2f Sv   " % (np.nanmin(sfm[i]), np.nanmax(sfm[i]))
        sys.stdout.write(u"\u001b[1000D" + "processing movie [%3d/%3d]" % (i+1,nT))
        sys.stdout.flush()
        ax.set_title(ptitle, fontsize = 12, y = 1.02)
        #
        return c, cf

    # https://stackoverflow.com/questions/18797175/animation-with-pcolormesh-routine-in-matplotlib-how-do-i-initialize-the-data
    # ani = animation.FuncAnimation(fig,update_img,  fargs = (nT,ax,im,c1,sfm,titlezer,) ,frames = nT, blit=False,repeat=False)
    anim = animation.FuncAnimation(fig, animate, frames=nT, blit=False, repeat=False)
    writer = animation.writers['ffmpeg'](fps=4)
    anim.save('%s.mp4' % (psave), writer=writer, dpi=200)
    plt.close("all")
    print("\nsaving : %s" % psave)
else :
    if save :
        print("\nsaving : %s" % psave)
        fig.savefig(psave)
    plt.show()
##################################################
# plt.show()
