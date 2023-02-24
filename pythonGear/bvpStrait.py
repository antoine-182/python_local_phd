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
parser = argparse.ArgumentParser(description = "Plot the porosity field across the strait in the ridge configuration")
parser.add_argument("netfile", type=str,
                    help="Netcdf file path")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
parser.add_argument("-s","--save", action="count", default=0,
                    help="save the current figure")

parser.add_argument("-m","--meshmask", type=str,
                    default="mesh_mask.nc",
                    help="meshmask associated to the file")

parser.add_argument("-dy", type=float, default = 5,
                    help="phi spacegrid")
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
rn_dy = args.dy # km
########################################################

# pmm = "ridge/zco/mesh_mask.nc"
# pdt = "ridge/zco/RIDGE_ref_zco_12h_grid_T.nc"
#
# psave = "zer.png"
# save = 0
# rn_dy = 5. # km

#
# rn_hdepth   = 5500.   ; rn_hsill   = 4430. # m  m
# rn_height   = 3180.   ; rn_width   = 5     # m  km
# rn_stddev   = 200.    ; rn_smt     = 5.    # km km
#
rn_hdepth   = 5500.   ; rn_hsill   = 4500. # m  m
rn_height   = 3000.   ; rn_width   = 30    # m  km
rn_stddev   = 200.    ; rn_smt     = 5     # km km
rn_rstd = 6.
#zphi0 = 1000.*rn_dy/2.
zphi0=0.
#

""" *****************************************************************
"""
# thetao = dtT.variables['toce']
psave = "strait.png"

""" Mask
"""
mm = nc4.Dataset(pmm)
tmask = mm.variables['tmask'][0,:,:,:]

nK,nJ,nI = np.shape(tmask)
midY = nJ//2
midX = nI//2
midZ = nK//2

tmask = tmask[:,:,midX]
gphit = mm.variables['gphit'][0,:,midX]
gphiv = mm.variables['gphiv'][0,:,midX]

""" Porosity
"""
dtT = nc4.Dataset(pdt)
try:
    rpoT = dtT.variables['rpot'][0,:,1,:]
except:
    rpoT = tmask

""" Geometrical setting
"""

# depth
ztmask =             tmask [  :,midY     ]
depw = mm.variables['gdepw_0'][0,:,midY,midX]
dept = mm.variables['gdept_0'][0,:,midY,midX]
# hdepth = np.max(np.cumsum(e3t*depthmask))
# print(hdepth)

# e3w = mm.variables['e3w_0'][0,:,:,midX]
e3t = mm.variables['e3t_0'][0,:,:,midX]
e2t = mm.variables['e2t'  ][0,  :,midX]

# gridk et gridy sont les sommets définissant les coins des pixels de pcolormesh
# e3w[0,:] /=2. # premier niveau dépasse la surface de la mer
gridk = np.zeros((e3t.shape[0]+1, e3t.shape[1]+1))
gridk[1:,:-1] = np.cumsum(e3t,axis=0) # somme sur la vertical de i,j=0,0 (upper left corner)
                                      # vers les i ascendant (z descendant)
gridk[1:,-1]  = gridk[1:,-2]          # duplique la dernière colonne avec l'avant dernière

gridy = np.zeros(gridk.shape)               # grille y
gridy[:,:-1] = gphit    - e2t    /2         # sur toute la colonne,
gridy[:, -1] = gphit[-1]+ e2t[-1]/2

gridy = gridy/1E3 ; gridk=gridk/100
nx = np.shape(gridy)[0] ; nk = np.shape(gridy)[1]

##### there is no equivalent for glam in depths
gridy = np.zeros((nJ,nK))  ;  gridz = gdepth*0.
for i in range(1,nI):
    gridx[:,i] = glamu[:  ,i-1]
for j in range(1,nJ):
    gridy[j,:] = gphiv[j-1,:  ]
gridx[:,0] = glamu[:,0] - e1t[:,0] ; gridy[0,:] = gphiv[0,:] - e2t[0,:]


""" Map
    A good manner to test if the grid match the correct data, is to
    plot without masking the datafield.
    vmin and vmax are necessary to turn blank 0 porosity
"""
palette = plt.get_cmap('Blues')
# fig, ax = plt.subplots(figsize = [12, 8])
fig, ax = plt.subplots(dpi=200)

rpoT = np.ma.masked_where(tmask==0,rpoT)

Y = np.cumsum(e2t[:-1]/1E3) - .5 # ??

z1d = 0.5 * ( np.tanh( (gphiv -zphi0 - 1E3*rn_width/2.) / (rn_smt*1E3) )  \
            - np.tanh( (gphiv -zphi0 + 1E3*rn_width/2.) / (rn_smt*1E3) )  )
zht =  rn_hdepth - ( rn_hdepth - ( rn_height + ( rn_height - rn_hsill ) * z1d ) )
zht = zht/100

im = plt.pcolormesh(gridy, gridk, rpoT,
                   vmin=0,vmax=1,
                   cmap = palette)
ax.plot(gphit/1E3,zht, color='tab:red', marker='o',
        linewidth = 0.8, markersize = 1.5)

# cbar = plt.colorbar(im)

ax.plot(gridy.T,gridk.T, 'w-', lw=0.3)
ax.plot(gridy  ,gridk  , 'w-', lw=0.3)

ax.patch.set_color('peru')
ax.set_ylabel("Z (100m)")
ax.set_xlabel("Y (km)")

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", width=1.5, length = 7, labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', width=1.5, length = 7, labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)

plt.tight_layout()
# ax.set_aspect(aspect=1) # data coordinate 'equal'

##################################################


x1=-4*rn_dy;x2=+4*rn_dy;y1=20;y2=50
# # x1=35.5;x2=42;y1=44;y2=52.5
ax.set_ylim(y2,y1)
ax.set_xlim(x1,x2)
#
# for i in range(int(x1),int(x2)):
#     for j in range(int(y1),int(y2)):
#         if np.ma.is_masked(rpoT[j,i]):
#             continue
#         ax.text(gridk[i,j]+0.5,gridy[i,j]+0.5,
#                 r"%2.0f" % (rpoT[j,i]*100)+"%",
#                 ha="center", va="center", fontsize = 6,color='black' )

if save :
    print("\nsaving : %s" % psave)
    fig.savefig(psave)
    plt.close()

plt.show()

##################################################
