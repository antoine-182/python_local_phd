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

import argparse
parser = argparse.ArgumentParser(description = "Plot the porosity field in the overflow configuration")
parser.add_argument("netfile", type=str,
                    help="Netcdf file path")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
parser.add_argument("-s","--save", action="count", default=0,
                    help="save the current figure")

parser.add_argument("-m","--meshmask", type=str,
                    default="mesh_mask.nc",
                    help="meshmask associated to the file")

args = parser.parse_args()

if args.verbose:
    print("Running '{}'".format(__file__))
    print("netcdf used " + args.netfile)


pdt = args.netfile
try:
    pmm = args.meshmask
    nc4.Dataset(pmm)
except:
    try:
        pmm = "/Users/gm/Documents/pythonGear/meshmaskOVERFLOW/"+args.meshmask
        nc4.Dataset(pmm)
    except:
        exit
save = args.save
psave = "rpoover.png"
########################################################

# pmm = "meshmaskoverflow/mesh_mask_x2z2.nc"
# pmm = "meshmaskoverflow/mesh_mask.nc"
# # pdt = "OVF_bvp_x2z2_grid_T.nc"
# pdt = "OVF_bvp_x2z2_grid_T.nc"
#
# psave = "zer.png"
# save = 0

""" *****************************************************************
"""
dtT = nc4.Dataset(pdt)
rpoT = dtT.variables['rpot'][0,:,1,:]
""" Mask
"""
mm = nc4.Dataset(pmm)
tmask = mm.variables['tmask'][0,:,1,:]

nI,nJ = np.shape(tmask)

""" Geometrical setting
"""
e3w = mm.variables['e3w_0'][0,:,1,:]
e1u = mm.variables['e1u'][0,1,:]

# gridk et gridx sont les sommets définissant les pixels de pcolormesh
gridk = np.zeros((e3w.shape[0]+1, e3w.shape[1]+1))
gridk[1:,:-1] = np.cumsum(e3w,axis=0) # somme sur la vertical de i,j=0,0 (upper left corner)
                                      # vers les i ascendant (z descendant)
gridk[1:,-1]  = gridk[1:,-2]          # duplique la dernière colonne avec l'avant dernière

gridx = np.zeros((e3w.shape[0]+1, e3w.shape[1]+1))
gridx[:,1:] = np.repeat([np.cumsum(e1u)], repeats = gridk.shape[0],axis = 0)

gridx = gridx/1E3 ; gridk = gridk/20.
nx = np.shape(gridx)[0] ; nk = np.shape(gridx)[1]


""" Map
    A good manner to test if the grid match the correct data, is to
    plot without masking the datafield.
    vmin and vmax are necessary to turn blank 0 porosity
"""
palette = plt.get_cmap('Blues')
# fig, ax = plt.subplots(figsize = [12, 8])
fig, ax = plt.subplots(dpi=200)

thetao = dtT.variables['thetao_inst']
rpoT = np.ma.masked_where(thetao[0,:,1,:]<5.,rpoT)

X = np.cumsum(e1u[:-1]/1E3)-.5

zht = 500. + 750.*( 1 + np.tanh( (X-40)/7 ) )
zht = zht / 20.
im = plt.pcolormesh(gridx, gridk, rpoT,
                   vmin=0,vmax=1,
                   cmap = palette)
ax.plot(X+1.,zht, color='tab:red', marker='o',
        linewidth = 0.8, markersize = 1.5)

# cbar = plt.colorbar(im)

# ax.plot(gridx  ,gridy  , 'w-', lw=0.4, zorder = 1)
# ax.plot(gridx.T,gridy.T, 'w-', lw=0.4, zorder = 1)
# ax.plot(glamt  ,gphit  , 'w--', lw=0.3, zorder = 1)
# ax.plot(glamt.T,gphit.T, 'w--', lw=0.3, zorder = 1)
ax.plot(gridx.T,gridk.T, 'w-', lw=0.3)
ax.plot(gridx  ,gridk  , 'w-', lw=0.3)

ax.patch.set_color('peru')
ax.set_ylabel("k-level")
ax.set_xlabel("length (km)")

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", width=1.5, length = 7, labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', width=1.5, length = 7, labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)

plt.tight_layout()
ax.set_aspect(aspect=1) # data coordinate 'equal'

##################################################


# x1=51;x2=65;y1=97;y2=101
x1=35.5;x2=42;y1=44;y2=52.5
ax.set_ylim(y2,y1)
ax.set_xlim(x1,x2)

for i in range(int(x1),int(x2)):
    for j in range(int(y1),int(y2)):
        if np.ma.is_masked(rpoT[j,i]):
            continue
        ax.text(gridk[i,j]+0.5,gridx[i,j]+0.5,
                r"%2.0f" % (rpoT[j,i]*100)+"%",
                ha="center", va="center", fontsize = 6,color='black' )

if save :
    print("\nsaving : %s" % psave)
    fig.savefig(psave)

##################################################
