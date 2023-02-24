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


# https://stackoverflow.com/questions/3940128/how-can-i-reverse-a-list-in-python
def Atick(A,B,delta,M=0):
    # A ... M ... B every dx
    #        List from A to M, reversed                         List from M to B and M is
    #                                                           taken away as it already exist
    return([ 2*M-x for x in np.arange(M,np.abs(2*M-A),delta)[::-1] ] + [x for x in np.arange(M,B+delta,delta)[1:] ])

""" *****************************************************************
"""

import argparse
parser = argparse.ArgumentParser(description = "Plot the porosity field for AM98")
parser.add_argument("netfile", type=str,
                    help="Netcdf files path (must give a pattern)")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
parser.add_argument("-d", type=int,
                    default=0,
                    help="lenght (number of cells -gc) of the penalisation")
parser.add_argument("-n", type=int,
                    default=4,
                    help="grid spacing in km")
parser.add_argument("-dx", type=real,
                    default=25*1000,
                    help="grid spacing in km (default 25km)")
parser.add_argument("-t","--theta", type=real,
                    default=0.,
                    help="rotation angle of the grid (default 0°)")
parser.add_argument("-s","--save", action="count", default=0,
                    help="save the current figure")

parser.add_argument("-m","--meshmask", type=str,
                    default="mesh_mask_025_0.nc",
                    help="meshmask associated to the file")

args = parser.parse_args()

if args.verbose:
    print("Running '{}'".format(__file__))
    print("netcdf used " + args.netfile)

""" *****************************************************************
"""
theta = float(args.theta)
dx = 100E3 # 100km
dx/= float(args.n)
dy=dx
print(args.meshmask)
pdtT,pdtU,pdtV,pdtF = [args.netfile + "%s.nc" % (s) for s in ["T","U","V","F"]]
psave = "rpo.png"
try:
    pmm = args.meshmask
    nc4.Dataset(pmm)
except:
    try:
        pmm = "/Users/gm/Documents/pythonGear/meshmaskAM98/"+args.meshmask
        nc4.Dataset(pmm)
    except:
        exit
""" *****************************************************************
"""


""" *****************************************************************
"""

# ssh
dtT = nc4.Dataset(pdtT)
dtU = nc4.Dataset(pdtU)
dtV = nc4.Dataset(pdtV)
dtF = nc4.Dataset(pdtF)
rpoT = dtT.variables['rpo'][-1,:,:]
rpoU = dtU.variables['rpou'][-1,:,:]
rpoV = dtV.variables['rpov'][-1,:,:]

try:
    rpoF = dtT.variables['rpof'][-1,:,:]
except:
    try:
        rpoF = dtF.variables['rpof'][-1,:,:]
    except:
        rpoF = np.nan*np.copy(rpoT)
""" Mask
"""
mm = nc4.Dataset(pmm)
tmask = mm.variables['tmask'][0,0]
umask = mm.variables['umask'][0,0]
vmask = mm.variables['vmask'][0,0]

nI,nJ = np.shape(tmask)
ssfmask = 0.*np.copy(tmask)
for i in range(nI-1):
    for j in range(nJ-1):
        ssfmask[i,j] = np.nanmax([tmask[i  ,j+1], tmask[i+1,j+1],   \
                                  tmask[i  ,j  ], tmask[i+1,j  ]    ])

""" Geometrical setting
"""
# glamt, gphit
glamt  = dtT.variables['nav_lon'][:] # x
gphit  = dtT.variables['nav_lat'][:] # y

glamu  = dtU.variables['nav_lon'][:]
gphiu  = dtU.variables['nav_lat'][:]

glamv  = dtV.variables['nav_lon'][:]
gphiv  = dtV.variables['nav_lat'][:]

# grid centré sur les corners bas gauche
lx = dx*np.cos((45+theta)*np.pi/180)/np.sqrt(2)
ly = dx*np.sin((45+theta)*np.pi/180)/np.sqrt(2)
# ballec cotes droit/haut - ils seront cropés
nx = np.shape(glamt)[0] ; ny = np.shape(glamt)[1]
gridx = np.zeros((nx,ny))
gridx = glamt - lx
gridy = np.zeros((nx,ny))
gridy = gphit - ly

""" Map
    A good manner to test if the grid match the correct data, is to
    plot without masking the datafield.
    vmin and vmax are necessary to turn blank 0 porosity
"""
palette = plt.get_cmap('Blues')
fig, ax = plt.subplots(dpi=200)

rpoT[tmask==0] = 0
# rpoT = smoother(rpoT, n = 0)

""" TO BE CHANGED
"""
#print("Warning : rpoU, rpoV field not read !!")
#rpoU = np.copy(rpoT) ; rpoV = np.copy(rpoT)
#rpoU = np.ma.masked_where(umask==0,rpoU)
#rpoV = np.ma.masked_where(vmask==0,rpoV)

rpoT = np.ma.masked_where(tmask==0,rpoT)
rpoU = np.ma.masked_where(umask==0,rpoU)
rpoV = np.ma.masked_where(vmask==0,rpoV)

try:
    rpoF = np.ma.masked_where(ssfmask==0,rpoF)
except:
    print('no rpoF')

im = ax.pcolormesh(gridx,gridy, rpoT,
                   vmin=0,vmax=1,
                   cmap = palette)
# cbar = plt.colorbar(im)

ax.plot(gridx  ,gridy  , 'w-', lw=0.4, zorder = 1)
ax.plot(gridx.T,gridy.T, 'w-', lw=0.4, zorder = 1)
ax.plot(glamt  ,gphit  , 'w--', lw=0.3, zorder = 1)
ax.plot(glamt.T,gphit.T, 'w--', lw=0.3, zorder = 1)

ax.patch.set_color('peru')
# a1 = -30000; a2 = 100000; a3 = 950000; a4 = 1050000
float(args.d)
a1 = -( 1.5 + float(args.d) )*dx ; a2 = (1.5 + float(args.d) )*dx
a3 = 1000000 - dx ; a4 = 1000000 + dx
print("box x: %03d km --- %04d km" % (a1/1E3, a2/1E3))
print("box y: %03d km --- %04d km" % (a3/1E3, a4/1E3))
ax.set_xlim(a1,a2)
ax.set_ylim(a3,a4)
# ax.vlines(x=0, ymin = a3, ymax=a4, color = "greenyellow", linewidth = 2)
##################################################
Xtick = Atick(a1,a2,dx)
ax.set_xticks(Xtick)
ax.set_xticklabels(["%.0f" % (x/1000) for x in Xtick ])
# ax.set_xticklabels(["-25","0","25","50","75","100"])
ax.set_xlabel("X (km)")

Ytick = Atick(a3,a4,dx, M = 1000000)
ax.set_yticks(Ytick)
ax.set_yticklabels(["%.0f" % (x/1000) for x in Ytick ])
ax.set_ylabel("Y (km)")
###############################################################
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", width=1.5, length = 7, labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', width=1.5, length = 7, labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
##################################################

# https://matplotlib.org/3.1.1/gallery/images_contours_and_fields/image_annotated_heatmap.html
textcolors=["black", "white"] ; ft = 4 # 6 par défualt
# Normalize the threshold to the images color range.
threshold = 0.9
if threshold is not None:
    threshold = im.norm(threshold)
else:
    threshold = im.norm(np.nanmax(rpoT))/2.

from  matplotlib.patches import Ellipse
for i in range(nx):
    for j in range(ny):
        if ( glamt[i,j] < a1*0.9 or glamt[i,j] > a2*1.1 or gphit[i,j] < a3 or gphit[i,j] > a4 ) :
            continue
        ax.text(glamt[i,j],gphit[i,j], r"$\phi_t$=%1.2f" % (0 if np.ma.is_masked(rpoT[i,j]) else rpoT[i,j]),
        #ax.text(glamt[i,j],gphit[i,j], r"$\phi_t$=%1.2f" % (rpoT[i,j]),     # perd l'intéret du mask
        # ax.text(glamt[i,j],gphit[i,j], r"$\phi_t$=%1.2f" % (tmask[i,j]),
                  ha="center", va="center",
                  fontsize = ft,color=textcolors[int(im.norm(0 if np.ma.is_masked(rpoT[i,j]) else rpoT[i,j]) > threshold)] )
#                  fontsize = ft,color='black' )

        # https://stackoverflow.com/questions/52056475/python-plot-rectangles-of-known-size-at-scatter-points
        height =  dx/10.                # m (hauteur)
for i in range(nx):
    for j in range(ny):
        if ( glamv[i,j] < a1*0.9 or glamv[i,j] > a2*1.1 or gphiv[i,j] < a3 or gphiv[i,j] > a4 ) :
            continue
        # multiplied by nan so the rectangle is erased
        widthv  = dx * (4./5.) * (np.nan if np.ma.is_masked(rpoV[i,j]) else rpoV[i,j]) # m (longueur)
        # la rotation de Rectangle() se fait depuis xy (après la translation de xy donc)
        # pour les rectangles xy est le lower-left corner
        wv =  widthv * np.cos(theta*np.pi/180) - height * np.sin(theta*np.pi/180)
        hv =  widthv * np.sin(theta*np.pi/180) + height * np.cos(theta*np.pi/180)
        # V facets
        ax.add_patch(Rectangle(xy=(glamv[i,j]-wv/2, gphiv[i,j]-hv/2),
                               width= widthv, height=height, angle = theta,
                               linewidth=1, color='red', fill=True, zorder=2))
        ax.text(glamv[i,j],gphiv[i,j], r"$\phi_v$=%1.2f" % (0 if np.ma.is_masked(rpoV[i,j]) else rpoV[i,j]),
                  ha="center", va="center", rotation = theta,
                  fontsize = ft,color="w" )

for i in range(nx):
    for j in range(ny):
        if ( glamu[i,j] < a1*0.9 or glamu[i,j] > a2*1.1 or gphiu[i,j] < a3 or gphiu[i,j] > a4 ) :
            continue
        # multiplied by nan so the rectangle is erased
        widthu  = dx * (4./5.) * (np.nan if np.ma.is_masked(rpoU[i,j]) else rpoU[i,j]) # m (longueur)
        # cos(90+theta) et sin(90+theta)
        wu =   widthu * np.sin(theta*np.pi/180) + height * np.cos(theta*np.pi/180)
        hu = - widthu * np.cos(theta*np.pi/180) + height * np.sin(theta*np.pi/180)
        # U facet1
        ax.add_patch(Rectangle(xy=(glamu[i,j]-wu/2, gphiu[i,j]-hu/2),
                               width= widthu, height=height, angle = theta - 90,
                               linewidth=1, color='red', fill=True, zorder=2))
        ax.text(glamu[i,j],gphiu[i,j], r"$\phi_u$=%1.2f" % (0 if np.ma.is_masked(rpoU[i,j]) else rpoU[i,j]),
                  ha="center", va="center", rotation = theta - 90,
                  fontsize = ft,color="w" )

for i in range(nx):
    for j in range(ny):
        if ( glamt[i,j]+lx < a1*0.9 or glamt[i,j]+lx > a2*1.1 or gphit[i,j]+ly < a3 or gphit[i,j]+ly > a4 ) :
            continue
        # F point
        # https://stackoverflow.com/questions/52056475/python-plot-rectangles-of-known-size-at-scatter-points
        radiusf =  dx/10.             # m (hauteur)
        # multiplied by nan so the rectangle is erased
        radiusf  = dx * (1./8.) * (np.nan if (np.ma.is_masked(rpoF[i,j]) or rpoF[i,j]==0) else 1. ) # m (longueur)
        # pour les ellipse, xy est le centre de l'ellipse
        # V facets
        ax.add_patch(Circle(xy=(glamt[i,j] + lx, gphit[i,j] + ly),
                               radius = radiusf,
                               edgecolor = "black",
                               linewidth=1, facecolor='red', fill=True, zorder=3))
        ax.text(glamt[i,j] + lx,gphit[i,j] + ly, r"%1.2f" % (0 if np.ma.is_masked(rpoF[i,j]) else rpoF[i,j]),
        # ax.text(glamt[i,j] + lx,gphit[i,j] + ly, r"%1.2f" % (rpoF[i,j]),
                  ha="center", va="center",
                  fontsize = ft,color="w" )


# for i in range(nx):
#     for j in range(ny):
#         if ( gridx[i,j] < a1*0.9 or gridx[i,j] > a2*1.1 or gridy[i,j] < a3*0.9 or gridy[i,j] > a4*1.1 ) :
#             continue
#         ax.text(glamt[i,j],gphit[i,j], r"%.0f%%" % (0 if np.ma.is_masked(rpoT[i,j]) else (100*rpoT[i,j]) ),
#                   ha="center", va="center", rotation = 0,
#                   fontsize = 14,color=textcolors[int(im.norm(0 if np.ma.is_masked(rpoT[i,j]) else rpoT[i,j]) > threshold)] )
# #                  fontsize = ft,color='black' )

# https://matplotlib.org/examples/pylab_examples/dolphin.html
# circle = Circle((0, 0), 1, facecolor='none',
#                 edgecolor=(0, 0.8, 0.8), linewidth=3, alpha=0.5)
# ax.add_patch(circle)
# im.set_clip_path(circle)

# https://jdhao.github.io/2017/06/03/change-aspect-ratio-in-mpl/
ax.set_aspect(aspect=1) # data coordinate 'equal'
# https://stackoverflow.com/questions/9295026/matplotlib-plots-removing-axis-legends-and-white-spaces
# plt.axis('off')

##################################################

# plt.show()
# showsave(psave, fig, save)

# showsave(psave, fig, save)
if args.save :
    print("saving : %s" % psave)
    fig.savefig(psave, dpi = 200)
else :
    plt.show()
# fig.savefig(psave, facecolor=fig.get_facecolor())
# plt.close(fig)
