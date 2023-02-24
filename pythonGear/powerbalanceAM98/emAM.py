# -*- coding: utf-8 -*-

import netCDF4 as nc4
import numpy as np
import time, sys, os

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, FixedLocator, FixedFormatter,
                               NullLocator)
import matplotlib.gridspec as gridspec

import matplotlib.animation as animation
from pylab import *

""" *****************************************************************
"""

# import argparse
# parser = argparse.ArgumentParser(description = "Plot the energy balance at the end of the simulation for AM98")
# parser.add_argument("netfileU", type=str,
#                     help="Netcdf file path")
# parser.add_argument("netfileV", type=str,
#                     help="Netcdf file path")
# parser.add_argument("-v", "--verbose", action="store_true",
#                     help="increase output verbosity")
#
# parser.add_argument("-dx", type=real,
#                     default=25,
#                     help="grid spacing in km")
# parser.add_argument("-t","--theta", type=real,
#                     default=0.,
#                     help="rotation angle of the grid (default 0°)")
# parser.add_argument("-s","--save", action="count", default=0,
#                     help="save the current figure")
# parser.add_argument("-f","--film", action="count", default=0,
#                     help="make a film")
#
# parser.add_argument("-min","--minimum", type=float, default = 200,
#                     help="minimal value of the colorbar")
# parser.add_argument("-max","--maximum", type=float, default = 800,
#                     help="maximal value of the colorbar")
#
# parser.add_argument("-m","--meshmask", type=str,
#                     default="mesh_mask_025_0.nc",
#                     help="meshmask associated to the file")
#
# args = parser.parse_args()
#
# if args.verbose:
#     print("Running '{}'".format(__file__))
#     print("netcdf used " + args.netfile)

""" *****************************************************************
"""
# theta = float(args.theta)
# dx =  float(args.dx) ; dy =  float(args.dx) # 25 km
# dx*=1000;dy*=1000
# print(args.netfileU) ; print(args.netfileV)
# pdu = args.netfileU ; pdv = args.netfileV
#
# print(args.meshmask)
# # Result in, if there is no min/max declared, =None and pcolormesh handle this
# vmin = args.minimum ; vmax = args.maximum
#
# try:
#     pmm = args.meshmask
#     nc4.Dataset(pmm)
# except:
#     try:
#         pmm = "../"+args.meshmask
#         nc4.Dataset(pmm)
#     except:
#         try:
#             pmm = "/Users/gm/Documents/pythonGear/meshmaskAM98/"+args.meshmask
#             nc4.Dataset(pmm)
#         except:
#             exit

""" *****************************************************************
"""
# pdt = "../../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_REF_16/freeslip_0/AM98_ref_16_freeslip_0_5d_00010101_00051230_grid_T.nc"
# pdu = "../../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_REF_16/freeslip_0/AM98_ref_16_freeslip_0_10y_5d_00010101_00101230_grid_U.nc"
# pdv = "../../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_REF_16/freeslip_0/AM98_ref_16_freeslip_0_10y_5d_00010101_00101230_grid_V.nc"
# pmm = "../../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_REF_16/freeslip_0/mesh_mask.nc"
# dx = 6.25 ; dy =  6.25 # 25 km
# dx*=1000;dy*=1000
# theta = 0.
# vmin = -0.25 ; vmax = 0.25

pdt = "../../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_ref8/freeslip0/AM98_ref_8_freeslip_0_5d_00010101_00051230_grid_T.nc"
pdu = "../../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_ref8/freeslip0/AM98_ref_8_freeslip_0_5d_00010101_00051230_grid_U.nc"
pdv = "../../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_ref8/freeslip0/AM98_ref_8_freeslip_0_5d_00010101_00051230_grid_V.nc"
pmm = "../../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_ref8/freeslip0/mesh_mask.nc"
dx = 12.5 ; dy =  12.5 # 25 km
dx*=1000;dy*=1000
theta = 0.
# vmin = -0.25 ; vmax = 0.25

# pdt = "../../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_ref4/freeslip0/AM98_ref_4_freeslip_0_5d_00010101_00051230_grid_T.nc"
# pdu = "../../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_ref4/freeslip0/AM98_ref_4_freeslip_0_5d_00010101_00051230_grid_U.nc"
# pdv = "../../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_ref4/freeslip0/AM98_ref_4_freeslip_0_5d_00010101_00051230_grid_V.nc"
# pmm = "../../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_ref4/freeslip0/mesh_mask.nc"
# dx = 25 ; dy =  25 # 25 km
# dx*=1000;dy*=1000
# theta = 0.
# vmin = -0.25 ; vmax = 0.25
vmin = -1. ; vmax = 1.


dtt = nc4.Dataset(pdt)
dtu = nc4.Dataset(pdu)
dtv = nc4.Dataset(pdv)


""" grid construction
"""
# if args.film:
#     nT,_,_ = np.shape(ssh)
mm = nc4.Dataset(pmm)
tmask = mm.variables['tmask'][0,0]
glamt = mm.variables['glamt'][0,:,:]
gphit = mm.variables['gphit'][0,:,:]

nI,nJ = np.shape(tmask)

# grid centré sur les corners bas gauche
lx = dx*np.cos((45+theta)*np.pi/180)/np.sqrt(2)
ly = dx*np.sin((45+theta)*np.pi/180)/np.sqrt(2)
# ballec cotes droit/haut - ils seront cropés
nx = np.shape(glamt)[0] ; ny = np.shape(glamt)[1]
gridx = np.zeros((nx,ny))
gridx = glamt - lx
gridy = np.zeros((nx,ny))
gridy = gphit - ly

limx = dx ; limy = dx

""" data construction
"""
rho0 = 1000 #kg/m3
rfr = 1e-7 #  s-1
tau = 0.2  #  N/m2 or
nu = 500 # m2/s

ht = dtt.variables["e3t"][-1,0,:,:]
hu = dtu.variables["hu"][-1,:,:]       ; hv = dtv.variables["hv"][-1,:,:]
uu = dtu.variables["ssu"][-1,:,:]      ; vv = dtv.variables["ssv"][-1,:,:]
tu = dtu.variables["sozotaux"][-1,:,:] ; tv = dtv.variables["sometauy"][-1,:,:]
divgrad = np.zeros((nI,nJ))
""" Conventional """
# dxu = np.zeros((nI,nJ)) ; dxv = np.zeros((nI,nJ))
# dyu = np.zeros((nI,nJ)) ; dyv = np.zeros((nI,nJ))
# hf = np.zeros((nI,nJ))  ; fmask = np.zeros((nI,nJ))
# for i in range(1,nI-1):
#     for j in range(1,nJ-1):
#         dxu[i,j] = (uu[i  ,j  ] - uu[i-1,j  ]) / dx  # T
#         dyu[i,j] = (uu[i  ,j+1] - uu[i  ,j  ]) / dx  # F
#         dxv[i,j] = (vv[i+1,j  ] - vv[i  ,j  ]) / dx  # F
#         dyv[i,j] = (vv[i  ,j  ] - vv[i  ,j-1]) / dx  # T
#         hf[i,j] = 0.25*( ht[i  ,j  ]*tmask[i  ,j  ] \
#                        + ht[i+1,j  ]*tmask[i+1,j  ] \
#                        + ht[i  ,j+1]*tmask[i  ,j+1] \
#                        + ht[i+1,j+1]*tmask[i+1,j+1])
#         fmask[i,j] = tmask[i,j]*tmask[i+1,j]*tmask[i,j+1]*tmask[i+1,j+1]
# for i in range(1,nI-1):
#     for j in range(1,nJ-1):
#         divgrad[i,j] = (tmask[i+1,j]*ht[i+1,j  ]*dxu[i+1,j  ] - tmask[i,j]    *ht[i  ,j  ]*dxu[i  ,j  ]           \
#                       + fmask[i,j]  *hf[i,j]    *dyu[i,j]     - fmask[i  ,j-1]*hf[i  ,j-1]*dyu[i  ,j-1])*uu[i,j]  \
#                       +(fmask[i,j]  *hf[i,j]    *dxv[i,j]     - fmask[i-1,j  ]*hf[i-1,j  ]*dxv[i-1,j  ]           \
#                       + tmask[i,j+1]*ht[i  ,j+1]*dyv[i  ,j+1] - tmask[i,j]    *ht[i  ,  j]*dyv[i  ,j  ])*vv[i,j]
# divgrad/=dx
""" divrot """
curl = np.zeros((nI,nJ)) ; dive = np.zeros((nI,nJ))
hf = np.zeros((nI,nJ))  ; fmask = np.zeros((nI,nJ))
for i in range(1,nI-1):
    for j in range(1,nJ-1):
        # hf[i,j] = 0.25*( ht[i  ,j  ]*tmask[i  ,j  ] \
        #                + ht[i+1,j  ]*tmask[i+1,j  ] \
        #                + ht[i  ,j+1]*tmask[i  ,j+1] \
        #                + ht[i+1,j+1]*tmask[i+1,j+1])   # F+
        hf[i,j] = 0.25*( ht[i  ,j  ] \
                       + ht[i+1,j  ] \
                       + ht[i  ,j+1] \
                       + ht[i+1,j+1])   # F+
        dive[i,j] = nu * (hu[i+1,j  ]*uu[i+1,j  ] - hu[i,j]*uu[i,j]) \
                       + (hv[i  ,j+1]*vv[i  ,j+1] - hv[i,j]*vv[i,j]) / ht[i,j]  # T
        curl[i,j] = nu * (uu[i  ,j+1] - uu[i,j]) \
                       - (vv[i+1,j  ] - vv[i,j]) * hf[i,j]                     # F+
        fmask[i,j] = tmask[i,j]*tmask[i+1,j]*tmask[i,j+1]*tmask[i+1,j+1]  # freeslip
curl = curl*fmask ; dive = dive*tmask
for i in range(1,nI-1):
    for j in range(1,nJ-1):
        ze1 = ( - (curl[i  ,j  ] - curl[i,j-1] )            \
                +  dive[j+1,j  ] - dive[i,j  ] ) *uu[i,j]    # U (acc x hu)
        ze2 = (   (curl[i  ,j  ] - curl[i-1,j] )            \
                +  dive[j  ,j+1] - dive[i  ,j] ) *vv[i,j]   # V (acc x hv)
        divgrad[i,j] = (ze1 + ze2)/(dx*dx)
# divgrad=divgrad*tmask
# divgrad = tmask*divgrad
lap = rho0 * divgrad
""" linear friction """
wind = np.zeros((nI,nJ)) ; fric = np.zeros((nI,nJ)) ; lap = np.zeros((nI,nJ))
z2d = hu*uu*uu + hv*vv*vv
fric = - rho0 * rfr * z2d
""" wind """
z2d = tu*uu + tv*vv
wind = z2d
""" non linear advection """


# datap = np.ma.masked_where(tmask==0,datap)
# datas = np.ma.masked_where(tmask==0,datas)
data = wind + fric + lap
# print("power balance : %f W/m2" % (np.nansum(data)))
print("power balance : %.1f GW" % (np.nansum(data)*dx*dx/1E9))

""" data map (end experiment)
"""
psave = "em.png"
# palette = plt.get_cmap('RdYlBu_r')
palette = plt.get_cmap('seismic')
level = np.arange(vmin,vmax,0.02)

""" figure bordelum """
fig, ax = plt.subplots(dpi=200)

im = ax.pcolormesh(gridx, gridy, data,
                    vmin=vmin, vmax=vmax,
                    # levels = level,
                    cmap = palette)
# c0 = ax.contour(glamt,gphit, data,
#                     levels = level,
#                     vmin=vmin,vmax=vmax,
#                     linewidths =0.8, colors=('k',),linestyles = "solid")
# ax.clabel(c0, fmt='%.2f', colors='k', fontsize=10)

# im = ax.pcolormesh(gridx, gridy, datas,
#                     vmin=vmin, vmax=vmax, alpha = 0.8,
#                     # levels = level,
#                     cmap = palette)
# c1 = ax.contour(glamt,gphit, datas,
#                     levels = level,
#                     vmin=vmin,vmax=vmax,
#                     linewidths =0.8, colors=('b',),linestyles = "solid")
# ax.clabel(c1, fmt='%.2f', colors='k', fontsize=10)
#
# im = ax.pcolormesh(gridx, gridy, datap,
#                     vmin=vmin, vmax=vmax, alpha = 0.8,
#                     # levels = level,
#                     cmap = palette)
# c2 = ax.contour(glamt,gphit, datap,
#                     levels = level,
#                     vmin=vmin,vmax=vmax,
#                     linewidths =0.8, colors=('w',),linestyles = "dashed")
# ax.clabel(c2, fmt='%.2f', colors='k', fontsize=10)

cbar = plt.colorbar(im)
# if args.film==0 :
#     c1 = ax.contour(glamt,gphit, data,
#                     # levels = levelc,
#                     vmin=vmin,vmax=vmax,
#                     linewidths =0.3, colors=('k',),linestyles = "solid")
#     ax.clabel(c1, fmt='%3.0f', colors='k', fontsize=10)

cbar.set_label(r"Local Power balance - (W/m2)")


ax.set_xlim(-limx,2000000)
ax.set_xticks([0,500000,1000000,1500000,2000000])
ax.set_xticklabels(["0","500","1000","1500","2000"])
ax.set_xlabel("X (km)")

ax.set_ylim(-limy,2000000)
ax.set_yticks([0,500000,1000000,1500000,2000000])
ax.set_yticklabels(["0","500","1000","1500","2000"])
ax.set_ylabel("Y (km)")

titlezer  = "Power Balance = %+.1f GW \n" % (np.nansum(data)*dx*dx/1E9)
titlezer += "min = %3.1f m   max = %3.1f m   " % (np.min(data), np.max(data))
ax.set_title(titlezer,
              fontsize = 10, y = 1.02)
##################################################

ax.patch.set_color('0.')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", which = 'both', width=1.5, labelsize = 10, pad = 5, left = True, right = True)
ax.tick_params(axis = 'x', which = 'both', width=1.5, labelsize = 10, pad = 10, bottom = True, top = True)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
ax.set_aspect(aspect=1) # data coordinate 'equal'


# if args.film :
#     """ movie thing """
#     def update_img(n,nT,ax,im,theta,ptitle):
#         import sys
#         data = ssh[n,:,:]
#         data = np.ma.masked_where(tmask==0,data)
#         ptitle = '%s - %02d/%02d \n'%(pdtT,n+1,nT)
#         ptitle += "min = %3.1f m   max = %3.1f m   " % (np.min(data), np.max(data)) \
#                +  "$\Delta$ = %3.1f m" % (np.max(data) - np.min(data))
#         sys.stdout.write(u"\u001b[1000D" + "processing [%3d/%3d]" % (n+1,nT))
#         sys.stdout.flush()
#         ax.set_title(ptitle, fontsize = 10, y = 1.02)
#         # new_color = palette(data.T.ravel())
#         #
#         # for tp in c1.collections:
#         #     tp.remove()
#         # c1 = ax.contour(glamt,gphit, data,
#         #                 # levels = levelc,
#         #                 vmin=vmin,vmax=vmax,
#         #                 linewidths =0.3, colors=('k',),linestyles = "solid")
#         # ax.clabel(c1, fmt='%3.0f', colors='k', fontsize=10)
#         #
#         # im.update({'facecolors':new_color})
#         data = data[:-1, :-1]
#         # ravel() converts C to a 1d-array
#         im.set_array(data.ravel())
#         return
#
#     # https://stackoverflow.com/questions/18797175/animation-with-pcolormesh-routine-in-matplotlib-how-do-i-initialize-the-data
#     ani = animation.FuncAnimation(fig,update_img,  fargs = (nT,ax,im,ssh,titlezer,) ,frames = nT, blit=False,repeat=False)
#     writer = animation.writers['ffmpeg'](fps=nT/10)
#     ani.save(psave,writer=writer,dpi=200)
#     print("\nsaving : %s" % psave)
# else :
#     if args.save :
#         print("\nsaving : %s" % psave)
#         fig.savefig(psave)
#         plt.close()
#     plt.show()
# plt.show()
print("\nsaving : %s" % psave)
fig.savefig(psave)
plt.close()
