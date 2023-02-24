#!/usr/bin/env python
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

""" *****************************************************************
"""

import argparse
parser = argparse.ArgumentParser(description = "Plot the e3 weighted density histogram of ridge test case")
parser.add_argument('-n','--netfile', nargs='+',
                    type=str, required=True,
                    help="Netcdf file path")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
parser.add_argument("-eps","--epsilon", type=float,
                    default = 1E-2,
                    help="tolerance on temperature")
parser.add_argument("-dT", type=float,
                    default=0.05,
                    help="width of a bin")
parser.add_argument("-j","--yjump", type=float,
                    default=0.15,
                    help="jump on y axis")
parser.add_argument("-s","--save", action="count", default=0,
                    help="save the current figure")
parser.add_argument("-ps","--psave", type=str,
                    default="zer",
                    help="name of the saved figure")
# parser.add_argument("-m","--meshmask", type=str,
#                     default="mesh_mask.nc",
#                     help="meshmask associated to the file")
parser.add_argument("-e1u","--e1u", type=float,
                    default = 1E3,
                    help="lateral space grid (instead of calling the meshmask)")

args = parser.parse_args()

if args.verbose:
    print("Running '{}'".format(__file__))
    print("netcdf used " + args.netfile)

""" *****************************************************************
"""
listpdt = args.netfile
save=args.save
eps=args.epsilon
dT=args.dT
Yjump=args.yjump


psave="histrhox"
if args.psave=="zer":
    import re
    for ndt in range(len(listpdt)):
        pdt=listpdt[ndt]
        m = re.search('OVF_(.+?)_grid', pdt)
        if m:
            psave+="_"+m.group(1)
        else:
            psave+="%d" % ndt
psave+=".png"

print(psave)

# try:
#     pmm = args.meshmask
#     nc4.Dataset(pmm)
# except:
#     try:
#         pmm = "/Users/gm/Documents/pythonGear/meshmaskOVERFLOW/"+args.meshmask
#         nc4.Dataset(pmm)
#     except:
#         exit
""" *****************************************************************
"""
# pmm = "meshmaskoverflow/mesh_mask_x2z2.nc"
# pmm = "meshmaskoverflow/mesh_mask.nc"
# pdt = "OVF_bvp_x2z2_grid_T.nc"
# listpdt = ["OVF_bvp_x2z2_grid_T.nc"]
#
# psave = "zer.png"
# save = 0
# eps = pow(10,-4)
# dT = 0.5 # not too low please
# Yjump=0.2

""" COUPE
"""

def splitaxis(ax0,ax1,yjump = 0.3):

    """ ylim """
    ax1.set_ylim(0,yjump)
    ax0.set_ylim(bottom = yjump, top=1.)
    """ xaxis and ticks """
    ax0.spines['bottom'].set_visible(False)
    ax1.spines['top'].set_visible(False)

    ax0.xaxis.tick_top()
    ax0.xaxis.set_minor_locator(AutoMinorLocator())
    ax0.tick_params(axis='x',which='both',
                    bottom=False,labelbottom=False,
                    labeltop=False)

    ax1.xaxis.tick_bottom()
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.tick_params(axis='x',which='both', top=False,labeltop=False)

    """ grid """
    for zer in [ax0,ax1]:
        zer.grid(which='major', linestyle='-', linewidth='0.3', color='black')
        zer.grid(which='minor', linestyle=':', linewidth='0.3', color='black')
        zer.xaxis.set_minor_locator(AutoMinorLocator())
        zer.yaxis.set_minor_locator(AutoMinorLocator())
        zer.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
        zer.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
        zer.tick_params(which='minor',length = 4)
        zer.tick_params(which='major',length = 6)
    return(ax0,ax1)


""" *****************************************************************
    Plot du test overflow - zps : z partial cells
    34 : time step
    101 : depth (z)
    3 : y dimension (1 cell embraced by two boundary layers)
    202 : x dimension
"""

""" Distribution of temperature """
rhobins = np.arange(28.,28.8,dT)
# Tbins = np.append(Tbins,20.)
""" Front position """
# Xframe = np.array([40,60,100])*1000 # in m
Xframe = np.array([60,100])*1000 # in m
# Xframe = np.array([10])*1000 # in m

colors=["royalblue","orangered","forestgreen"]
# colors=["orangered","forestgreen"]
# colors=["royalblue","tab:orange","limegreen"]
markers=["o","x","v","^"]

print(listpdt)
Ndt = len(listpdt)
distriT  = np.zeros( (Ndt, len(Xframe),len(Tbins) ) )
timeframe= np.zeros( (Ndt, len(Xframe)            ) )
# dt =  xr.open_dataset(pdt)
for ndt in range(Ndt):
    pdt=listpdt[ndt]

    dt = nc4.Dataset(pdt)
    thetao = dt.variables['thetao_inst']
    nT,_,_,nX = np.shape(thetao)

    theta = np.ma.masked_outside(thetao, 9.5,20.01)
    e3t = dt.variables['e3t_inst']
    # e1u = mm.variables['e1u'][0,1,:]
    # gridx = np.cumsum(e1u)
    gridx=np.cumsum(np.ones(np.shape(theta)[-1])*1E3)

    """ find correponding time frame """
    for ix in range(len(Xframe)):
        il = np.argmin([np.abs(Xframe[ix] - gridx[i]) for i in range(nX)])
        t=0 ; flag = True
        while ( flag and t<(nT-1) ):
            test = np.any(theta[t,:,1,il]<20.)
            if test:
                flag=False
            else:
                t+=1
        timeframe[ndt,ix]=t

    """ Initial volume of cold water
    do not forget the intermediate volume between 10 and 20 degree
    to bring them back to 10° volume
    """
    Vcold = 0
    for the in range(len(Tbins)):
        # print("[%1.2f,%1.2f[" % ((Tbins[the]-dT/2.)-eps, (Tbins[the]+dT/2.)) )
        tmp = np.ma.masked_outside(thetao[0,:,1,:],
                                   (Tbins[the]-dT/2.)-eps, (Tbins[the]+dT/2.))
        # volume d'eau à the
        e3 = np.ma.masked_where(tmp.mask, e3t[0,:,1,:])
        e3 = np.nansum(np.ma.filled(e3,0.))
        tmp2 = e3
        tmp2 *= (20.-Tbins[the])/10.
        Vcold+=tmp2
    # print("%f" % Vcold)
    """ Distribution over the range of temperature """
    for t in range(len(Xframe)):
        for the in range(len(Tbins)):
            tmp = np.ma.masked_outside(thetao[timeframe[ndt,t],:,1,:],
                                       (Tbins[the]-dT/2.)-eps, (Tbins[the]+dT/2.))
            # volume d'eau à the
            e3 = np.ma.masked_where(tmp.mask, e3t[timeframe[ndt,t],:,1,:])
            e3 = np.nansum(np.ma.filled(e3,0.))
            tmp2 = e3/Vcold
            tmp2 *= (20.-Tbins[the])/10.
            distriT[ndt,t,the] = tmp2
        # print(r"%3.0fkm (%.1fh) : %1.2f" % (Xframe[t]/1E3, timeframe[ndt,t]*0.5, np.sum(distriT[ndt,t,:])) )

""" figure bordelum """
""" Split figures
"""
# fig, ax = plt.subplots(2,1,figsize=(8,4), dpi=200)
# # at each front
# for t in range(len(Xframe)):
#     # on the upper and lower panel
#     for zer in ax:
#         # as many dt files
#         for ndt in range(Ndt):
#             # plot the distribution
#             zer.plot(Tbins, distriT[ndt,t,:],
#                     label=r"%3.0fkm (%.1fh)" % (Xframe[t]/1E3, timeframe[ndt,t]*0.5),
#                     marker=markers[ndt], color=colors[t],
#                     alpha=1., linewidth=.7, markersize=5.)
#
# ax[0],ax[1] = splitaxis(ax[0],ax[1], yjump=Yjump)
#
# ax[0].set_xlim(10.,20.)
# ax[1].set_xlim(10.,20.)
# ax[1].set_xlabel("Temperature (°C)")
#
# ax[0].legend(loc=9, ncol=len(Xframe))
#
# fig.add_subplot(ax[0])
# fig.add_subplot(ax[1])


""" figure bordelum """
""" Single figure
"""
fig, ax = plt.subplots(figsize=(8,4), dpi=200)
# at each front
for t in range(len(Xframe)):
    # as many dt files
    for ndt in range(Ndt):
        # plot the distribution
        ax.plot(Tbins, distriT[ndt,t,:],
                label=r"%3.0fkm (%.1fh)" % (Xframe[t]/1E3, timeframe[ndt,t]*0.5),
                marker=markers[ndt], color=colors[t],
                alpha=1., linewidth=.7, markersize=5.)

ax.grid(which='major', linestyle='-', linewidth='0.3', color='black')
ax.grid(which='minor', linestyle=':', linewidth='0.3', color='black')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
ax.set_xlim(10.,20.)
ax.set_xlabel("Temperature (°C)")
ax.set_ylim(0,0.3)
ax.set_ylabel("Dense water corresponding volume ratio")
ax.legend(loc=9, ncol=len(Xframe))

fig.add_subplot(ax)

#########

titlezer=""
for ndt in range(Ndt):
    pdt=listpdt[ndt]
    titlezer += "%s\n" % pdt.split('/')[-1]
titlezer += r"Dense water distribution (%dh simulation)" % (nT*0.5)
plt.suptitle(titlezer)

fig.subplots_adjust(top = 0.8, bottom=0.15, hspace = 0.02)

if save :
    print("saving : %s" % psave)
    fig.savefig(psave)
else :
    plt.show()


# for i in range(1,10):
#     eps = pow(10,-i)
#     theta_init =np.ma.masked_outside(thetao[0,:,1,:], 9.5,20.-eps)
#     e3 = np.ma.masked_where(theta_init.mask, e3t[0,:,1,:])
#     Vcold=np.sum(e3)
#     print(Vcold,i)
# 10832.418 1
# 11091.113 2 <- let's get this one
# 11283.476 3
# 11783.453 4
# 11863.453 5
# 12145.281 6
# 339954.22 7
# 339954.22 8
# 339954.22 9
