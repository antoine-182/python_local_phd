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

# import argparse
# parser = argparse.ArgumentParser(description = "Plot the streamfunction of a denser water mass (zcoor) with the rige config")
# parser.add_argument("netfile", type=str,
#                     help="Netcdf file path - U point")
# parser.add_argument("-v", "--verbose", action="store_true",
#                     help="increase output verbosity")
# parser.add_argument("-s","--save", action="count", default=0,
#                     help="save the current figure")
#
# parser.add_argument("-nt", type=int,
#                     default=-1,
#                     help="timeframe used")
# parser.add_argument("-f","--film", action="count", default=0,
#                     help="make a film")
#
# parser.add_argument("-dv", type=float, default = 0.1,
#                     help="color interval in the colorbar")
# parser.add_argument("-min","--minimum", type=float, default = -1.5,
#                     help="minimal value of the colorbar")
# parser.add_argument("-max","--maximum", type=float, default = 1.5,
#                     help="maximal value of the colorbar")
#
# parser.add_argument("-t","--pdt", type=str,
#                     default="",
#                     help="netfile on T point")
# parser.add_argument("-m","--meshmask", type=str,
#                     default="mesh_mask.nc",
#                     help="meshmask associated to the file")
#
# args = parser.parse_args()
#
# if args.verbose:
#     print("Running '{}'".format(__file__))
#     print("netcdf used " + args.netfile)
# ########################################################
#
# nt = args.nt
# pdu = args.netfile
# pdt = args.pdt
#
# # Result in, if there is no min/max declared, =None and pcolormesh handle this
# vmin = args.minimum ; vmax = args.maximum
# try:
#     pmm = args.meshmask
#     nc4.Dataset(pmm)
# except:
#     try:
#         pmm = "/Users/gm/Documents/pythonGear/meshmaskOVERFLOW/"+args.meshmask
#         nc4.Dataset(pmm)
#     except:
#         exit
#
# if args.dv==-1.:
#     dv = (vmax-vmin)/100.
# else:
#     dv = args.dv



########################################################
nt=-1

# pdt = "ridge/zco/RIDGE_ref_zco_12h_grid_T.nc"
# pmm = "ridge/zco/mesh_mask.nc"
# pdu = "ridge/zco/RIDGE_ref_zco_12h_grid_U.nc"

# pdt = "ridge/sco/RIDGE_ref_sco_12h_grid_T.nc"
# pmm = "ridge/sco/mesh_mask.nc"
# pdu = "ridge/sco/RIDGE_ref_sco_12h_grid_U.nc"

dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/"
pdt = dir+"EXP00/RIDGE_sco_12h_grid_T.nc"
pmm = dir+"EXP00/mesh_mask.nc"
pdu = dir+"EXP00/RIDGE_sco_12h_grid_U.nc"

save = 0
vmin = -1.5 ; vmax = 1.5 ; dv = 0.1
Ncolor=100

Ncolor=int((vmax-vmin)/dv)

"""
    Plot du test overflow - zps : z partial cells
    34 : time step
    101 : depth (z)
    3 : y dimension (1 cell embraced by two boundary layers)
    202 : x dimension
"""
dtu = nc4.Dataset(pdu)
dt  = nc4.Dataset(pdt)

uu  = dtu.variables['u_vol'][-1,:,:,:]   # volumic transport
toce = dt.variables['toce'] [-1,:,:,:]
nK,nY,nX = np.shape(uu)
nT=1

mm = nc4.Dataset(pmm)
umask = mm.variables['umask']
tmask = mm.variables['tmask']

uum  = np.ma.masked_where(umask==0,uu  )

print("domain size is (x,y) %dx%d with %d k-levels" % (nX,nY,nK))
middle = int(np.shape(toce)[-2] / 2)
print("middle is %d" % (middle))
# middle = 3

tmask = mm.variables['tmask'][0,:,:,:]

#######################################################

glamu  = mm.variables['glamu'][0,middle,:]
gdepth = mm.variables['gdepw_0'][0,:,middle,:]
# (x,y)=np.meshgrid(glamt,gdepth[:,0])
x = gdepth*0. ; y = gdepth
for k in range(nK):
    x[k,:]=glamu

########################################################

# UW points surrounded by T points
                   #    time, x, y
e3t = mm.variables['e3t_0'][0,:,middle,:]
e1t   = mm.variables['e1t'  ][0,middle,:]
glamt  = mm.variables['glamt'][0,middle,:]
########################################################

# gridk et gridx sont les sommets définissant les pixels de pcolormesh
gridk = np.zeros((e3t.shape[0]+1, e3t.shape[1]+1))
gridk[1:,:-1] = np.cumsum(e3t,axis=0) # somme sur la vertical de i,j=0,0 (upper left corner)
                                      # vers les i ascendant (z descendant)
gridk[1:,-1]  = gridk[1:,-2]          # duplique la dernière colonne avec l'avant dernière

gridx = np.zeros(gridk.shape)       # grille x
# gridx[:,:-1] = linx - e1u/2         # sur toute la colonne,
# gridx[:, -1] = linx[-1]+e1u[-1]/2

gridx[:,:-1] = glamt  #  - e1t    /2         # sur toute la colonne,
gridx[:, -1] = glamt[-1]+ e1t[-1]/2

# nK,_ = np.shape(gridk)
gridkp = np.copy(gridk)
gridxp = np.copy(gridx)
for n in range(nK):
    if (n%2!=0):
        gridkp[n,:] = np.nan
# gridkp = np.ma.masked_where(thetao[0,:,middle,:]<1.,gridkp)

########################################################

sfm = np.zeros((nT,nK,nX))
""" streamfunction """
for t in range(nT):
    for i in range(nX):
        for k in range(nK):
            sfm[t,k,i] += - np.sum(uum[k:,:,i])/1E6

########################################################
# if args.film :
#     data = sfm[0,:,:]
#     titlezer  = '%s - %02d/%02d \n'%(pdt,0+1,nT)
#     psave = "sfridge.mp4"
# else :
#     data = sfm[nt,:,:]
#     titlezer  = '%s - %02d/%02d \n'%(pdt,nt+1 if nt!=-1 else nT,nT)
#     psave = "sfridge.png"

titlezer  = '%s\n'%(pdt)
data = sfm[nt,:,:]

# data = np.ma.masked_where(umask[0,:,middle,:]==0, data)

########################################################
palette = plt.get_cmap('RdBu_r',Ncolor)
""" figure bordelum """
fig, ax = plt.subplots(figsize = [12, 8])

# plt.pcolormesh(gridx, gridk, gauss,alpha=1.,zorder = 0.)
im = plt.pcolormesh(gridx, gridk, data,
                    vmin=vmin, vmax=vmax,
                    cmap = palette)
c1 = plt.contour(x, y, data,
                 vmin=vmin, vmax=vmax,
                linewidths =1, colors=('k',), linestyles = "solid")
ax.clabel(c1, fmt='%1.1f', colors='k', fontsize=12)

cbar = plt.colorbar(im)
# cbar.set_label(r"Streamfunction $\psi$ (Sv)")
ax.plot(gridxp.T,gridkp.T, 'k-', lw=0.5)
# https://stackoverflow.com/questions/31877353/overlay-an-image-segmentation-with-numpy-and-matplotlib
# plt.pcolormesh(gridx, gridk, sill ,alpha=0.9,zorder = 0.)
# plt.pcolormesh(gridx, gridk, gauss,alpha=0.5,zorder = 0.)

titlezer += "min = %2.2f Sv   max = %2.2f Sv   " % (np.min(data), np.max(data))
ax.set_title(titlezer, fontsize = 18, y = 1.02)

ax.set_ylim(gridk[-1,0],0)
# ax.set_yticks([0,500,1000,1500,2000,2500])
# ax.set_yticklabels(["0","500","1000","1500","2000","2500"])
ax.set_ylabel("Z (m)")

# Xtick = np.arange(gridx[0,0],gridx[0,-1],100)
Xticks=np.round(glamt.data/1E3)*1E3
Xticks = Xticks[::3]
ax.set_xlim(gridx[0,0],gridx[0,-1])
ax.set_xticks(Xticks)
ax.set_xticklabels([ "%d"%(x/1E3) for x in Xticks])
ax.set_xlabel("X (km)")

ax.patch.set_color('0.')
# ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
plt.tight_layout(rect=[0,0,1,0.95])
# plt.tight_layout()

##################################################

# if args.film :
#     """ movie thing """
#     def update_img(n,nT,ax,im,theta,ptitle):
#         import sys
#         data = sfm[n,:,:]
#         data = np.ma.masked_where(umask[0,:,middle,:]==0, data)
#         ptitle = '%s - %02d/%02d \n'%(pdt,n+1,nT)
#         ptitle += "min = %2.2f Sv   max = %2.2f Sv   " % (np.min(data), np.max(data))
#
#         sys.stdout.write(u"\u001b[1000D" + "processing [%3d/%3d]" % (n+1,nT))
#         sys.stdout.flush()
#         ax.set_title(ptitle, fontsize = 18, y = 1.02)
#         im.set_array(data.ravel())
#         return
#
#     # https://stackoverflow.com/questions/18797175/animation-with-pcolormesh-routine-in-matplotlib-how-do-i-initialize-the-data
#     ani = animation.FuncAnimation(fig,update_img,  fargs = (nT,ax,im,theta,titlezer,) ,frames = nT, blit=False,repeat=False)
#     writer = animation.writers['ffmpeg'](fps=6)
#     ani.save(psave,writer=writer,dpi=200)
#     print("\nsaving : %s" % psave)
# else :
#     if args.save :
#         print("\nsaving : %s" % psave)
#         fig.savefig(psave)
#     plt.show()
##################################################
plt.show()
