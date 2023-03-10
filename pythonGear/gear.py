bathyOverflow.py                                                                                    000644  000766  000024  00000011071 14024346240 014031  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         #!/usr/bin/env python
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
parser = argparse.ArgumentParser(description = "plot the bathymetrie profil associated to the meshmask (2D only) in the overflow config")
parser.add_argument('-m','--meshmask', nargs='+',
                    type=str, required=True,
                    help="Netcdf file path")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
# parser.add_argument("-eps","--epsilon", type=float,
#                     default = 1E-2,
#                     help="tolerance on temperature")
# parser.add_argument("-dT", type=float,
#                     default=0.5,
#                     help="width of a bin")
# parser.add_argument("-j","--yjump", type=float,
#                     default=0.15,
#                     help="jump on y axis")
parser.add_argument("-s","--save", action="count", default=0,
                    help="save the current figure")
parser.add_argument("-ps","--psave", type=str,
                    default="zer",
                    help="name of the saved figure")
# parser.add_argument("-m","--meshmask", type=str,
#                     default="mesh_mask.nc",
#                     help="meshmask associated to the file")
# parser.add_argument("-e1u","--e1u", type=float,
#                     default = 1E3,
#                     help="lateral space grid (instead of calling the meshmask)")

args = parser.parse_args()

if args.verbose:
    print("Running '{}'".format(__file__))
    print("netcdf used " + args.netfile)

""" *****************************************************************
"""
listpmm = args.meshmask
save=args.save

psave="zer.png"
print(psave)

""" *****************************************************************
"""
# listpmm=['mesh_mask_bvp.nc']
# save=0

""" *****************************************************************
    Plot du test overflow - zps : z partial cells
    34 : time step
    101 : depth (z)
    3 : y dimension (1 cell embraced by two boundary layers)
    202 : x dimension
"""

colors=["royalblue","orangered","forestgreen"]
# colors=["orangered","forestgreen"]
# colors=["royalblue","tab:orange","limegreen"]
markers=["o","x","v","^"]

print(listpmm)
Nmm = len(listpmm)
# dt =  xr.open_dataset(pdt)

fig, ax = plt.subplots(figsize=(5,4), dpi=200)
# at each front

for nmm in range(Nmm):
    pmm=listpmm[nmm]
    if pmm!="mesh_mask_x0z0.nc":
        continue
    mm = nc4.Dataset(pmm)

    e3t = mm.variables['e3t_0'][0,:,1,:]
    Nz,Nx = np.shape(e3t)
    e1u = mm.variables['e1u'][0,1,:]
    gridx = np.cumsum(e1u)/1E3
    tmask= mm.variables['tmask'][0,:,1,:]

    e3 = np.ma.masked_where(tmask==0, e3t)
    h0 = np.nansum(np.ma.filled(e3,0.),axis=0)

for nmm in range(Nmm):
    pmm=listpmm[nmm]
    if pmm=="mesh_mask_x0z0.nc":
        continue
    mm = nc4.Dataset(pmm)

    e3t = mm.variables['e3t_0'][0,:,1,:]
    Nz,Nx = np.shape(e3t)
    e1u = mm.variables['e1u'][0,1,:]
    gridx = np.cumsum(e1u)/1E3
    tmask= mm.variables['tmask'][0,:,1,:]

    e3 = np.ma.masked_where(tmask==0, e3t)
    h = np.nansum(np.ma.filled(e3,0.),axis=0)

    ax.plot(gridx, h-h0,
            label="%s" % pmm.split('/')[-1],
            linewidth=0.9, alpha=0.9)



ax.grid(which='major', linestyle='-', linewidth='0.3', color='black')
ax.grid(which='minor', linestyle=':', linewidth='0.3', color='black')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)

ax.set_xlim(20,60)
ax.set_xlabel("distance (km)")
ax.set_ylim(-8,8)
ax.set_ylabel(r"$\Delta$(m)")

ax.legend(loc=1, fontsize=8)

fig.add_subplot(ax)

# titlezer=""
# for nmm in range(Nmm):
#     pmm=listpmm[nmm]
#     titlezer += "%s\n" % pdt.split('/')[-1]
# titlezer += ""
plt.suptitle("Difference with non penalised water column depth")

fig.subplots_adjust(top = 0.9, bottom=0.15, hspace = 0.02)

if save :
    print("saving : %s" % psave)
    fig.savefig(psave)
else :
    plt.show()
                                                                                                                                                                                                                                                                                                                                                                                                                                                                       bvpAM.py                                                                                            000644  000766  000024  00000025620 13763771670 012232  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         # -*- coding: utf-8 -*-

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

parser.add_argument("-dx", type=real,
                    default=25*1000,
                    help="grid spacing in km (default 25km)")
parser.add_argument("-t","--theta", type=real,
                    default=0.,
                    help="rotation angle of the grid (default 0??)")
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
dx =  args.dx ; dy =  args.dx # 25 km
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

# grid centr?? sur les corners bas gauche
lx = dx*np.cos((45+theta)*np.pi/180)/np.sqrt(2)
ly = dx*np.sin((45+theta)*np.pi/180)/np.sqrt(2)
# ballec cotes droit/haut - ils seront crop??s
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
a1 = -1.5*dx ; a2 = 4*dx
a3 = 1000000 - 2*dx ; a4 = 1000000 + 2*dx
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
textcolors=["black", "white"]
# Normalize the threshold to the images color range.
threshold = 0.9
if threshold is not None:
    threshold = im.norm(threshold)
else:
    threshold = im.norm(np.nanmax(rpoT))/2.

from  matplotlib.patches import Ellipse
for i in range(nx):
    for j in range(ny):
        if ( gridx[i,j] < a1*0.9 or gridx[i,j] > a2*1.1 or gridy[i,j] < a3*0.9 or gridy[i,j] > a4*1.1 ) :
            continue
        ax.text(glamt[i,j],gphit[i,j], r"$\phi_t$=%1.2f" % (0 if np.ma.is_masked(rpoT[i,j]) else rpoT[i,j]),
        #ax.text(glamt[i,j],gphit[i,j], r"$\phi_t$=%1.2f" % (rpoT[i,j]),     # perd l'int??ret du mask
        # ax.text(glamt[i,j],gphit[i,j], r"$\phi_t$=%1.2f" % (tmask[i,j]),
                  ha="center", va="center",
                  fontsize = 6,color=textcolors[int(im.norm(0 if np.ma.is_masked(rpoT[i,j]) else rpoT[i,j]) > threshold)] )
#                  fontsize = 6,color='black' )

        # https://stackoverflow.com/questions/52056475/python-plot-rectangles-of-known-size-at-scatter-points
        height =  dx/10.                # m (hauteur)
        # multiplied by nan so the rectangle is erased
        widthu  = dx * (4./5.) * (np.nan if np.ma.is_masked(rpoU[i,j]) else rpoU[i,j]) # m (longueur)
        widthv  = dx * (4./5.) * (np.nan if np.ma.is_masked(rpoV[i,j]) else rpoV[i,j]) # m (longueur)
        # la rotation de Rectangle() se fait depuis xy (apr??s la translation de xy donc)
        # pour les rectangles xy est le lower-left corner
        wv =  widthv * np.cos(theta*np.pi/180) - height * np.sin(theta*np.pi/180)
        hv =  widthv * np.sin(theta*np.pi/180) + height * np.cos(theta*np.pi/180)
        # cos(90+theta) et sin(90+theta)
        wu =   widthu * np.sin(theta*np.pi/180) + height * np.cos(theta*np.pi/180)
        hu = - widthu * np.cos(theta*np.pi/180) + height * np.sin(theta*np.pi/180)
        # V facets
        ax.add_patch(Rectangle(xy=(glamv[i,j]-wv/2, gphiv[i,j]-hv/2),
                               width= widthv, height=height, angle = theta,
                               linewidth=1, color='red', fill=True, zorder=2))
        ax.text(glamv[i,j],gphiv[i,j], r"$\phi_v$=%1.2f" % (0 if np.ma.is_masked(rpoV[i,j]) else rpoV[i,j]),
                  ha="center", va="center", rotation = theta,
                  fontsize = 6,color="w" )
        # U facet1
        ax.add_patch(Rectangle(xy=(glamu[i,j]-wu/2, gphiu[i,j]-hu/2),
                               width= widthu, height=height, angle = theta - 90,
                               linewidth=1, color='red', fill=True, zorder=2))
        ax.text(glamu[i,j],gphiu[i,j], r"$\phi_u$=%1.2f" % (0 if np.ma.is_masked(rpoU[i,j]) else rpoU[i,j]),
                  ha="center", va="center", rotation = theta - 90,
                  fontsize = 6,color="w" )
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
                  fontsize = 6,color="w" )


# for i in range(nx):
#     for j in range(ny):
#         if ( gridx[i,j] < a1*0.9 or gridx[i,j] > a2*1.1 or gridy[i,j] < a3*0.9 or gridy[i,j] > a4*1.1 ) :
#             continue
#         ax.text(glamt[i,j],gphit[i,j], r"%.0f%%" % (0 if np.ma.is_masked(rpoT[i,j]) else (100*rpoT[i,j]) ),
#                   ha="center", va="center", rotation = 0,
#                   fontsize = 14,color=textcolors[int(im.norm(0 if np.ma.is_masked(rpoT[i,j]) else rpoT[i,j]) > threshold)] )
# #                  fontsize = 6,color='black' )

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
    fig.savefig(psave)
else :
    plt.show()
# fig.savefig(psave, facecolor=fig.get_facecolor())
# plt.close(fig)
                                                                                                                bvpOverflow.py                                                                                      000644  000766  000024  00000011272 14024346270 013517  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         # -*- coding: utf-8 -*-

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

# gridk et gridx sont les sommets d??finissant les pixels de pcolormesh
gridk = np.zeros((e3w.shape[0]+1, e3w.shape[1]+1))
gridk[1:,:-1] = np.cumsum(e3w,axis=0) # somme sur la vertical de i,j=0,0 (upper left corner)
                                      # vers les i ascendant (z descendant)
gridk[1:,-1]  = gridk[1:,-2]          # duplique la derni??re colonne avec l'avant derni??re

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
                                                                                                                                                                                                                                                                                                                                      bvpStrait.py                                                                                        000644  000766  000024  00000012751 14027173167 013173  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         # -*- coding: utf-8 -*-

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

""" *****************************************************************
"""
# thetao = dtT.variables['toce']
psave = "strait.png"

""" Mask
"""
mm = nc4.Dataset(pmm)
tmask = mm.variables['tmask'][0,:,:,:]

nK,nJ,nI = np.shape(tmask)
midY = nJ//2 + 1
midX = nI//2 + 1

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

# gridk et gridy sont les sommets d??finissant les pixels de pcolormesh
# e3w[0,:] /=2. # premier niveau d??passe la surface de la mer
gridk = np.zeros((e3t.shape[0]+1, e3t.shape[1]+1))
gridk[1:,:-1] = np.cumsum(e3t,axis=0) # somme sur la vertical de i,j=0,0 (upper left corner)
                                      # vers les i ascendant (z descendant)
gridk[1:,-1]  = gridk[1:,-2]          # duplique la derni??re colonne avec l'avant derni??re

gridy = np.zeros(gridk.shape)               # grille y
gridy[:,:-1] = gphit    - e2t    /2         # sur toute la colonne,
gridy[:, -1] = gphit[-1]+ e2t[-1]/2

gridy = gridy/1E3 ; gridk=gridk/100
nx = np.shape(gridy)[0] ; nk = np.shape(gridy)[1]


""" Map
    A good manner to test if the grid match the correct data, is to
    plot without masking the datafield.
    vmin and vmax are necessary to turn blank 0 porosity
"""
palette = plt.get_cmap('Blues')
# fig, ax = plt.subplots(figsize = [12, 8])
fig, ax = plt.subplots(dpi=200)

rpoT = np.ma.masked_where(tmask==0,rpoT)

Y = np.cumsum(e2t[:-1]/1E3)-.5
#
rn_hdepth   = 5500.   ; rn_hsill   = 4430.  # m  m
rn_height   = 3180.   ; rn_width   = 20     # m  km
rn_stddev   = 200.    ; rn_smt     = 100.   # km m
#
zphi0 = 1000.*rn_dy/2.
# zphi0=0.
z1d = 0.5 * ( np.tanh( (gphiv -zphi0 - 1E3*rn_width/2.) / rn_smt )  \
            - np.tanh( (gphiv -zphi0 + 1E3*rn_width/2.) / rn_smt )  )
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
                       diffplotoverflow.py                                                                                 000644  000766  000024  00000010312 13761420745 014600  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         # -*- coding: utf-8 -*-

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
from pylab import *
# import cmocean

from useful import *
save = 1

""" *****************************************************************
    Plot du test overflow - zps : z partial cells
    34 : time step
    101 : depth (z)
    3 : y dimension (1 cell embraced by two boundary layers)
    202 : x dimension
"""


# """ e3w_cellcentered """
# pdtA = "e3w_cellcentered/overflow_zps_e3w_onoldTpoint_grid_T.nc"
# pmm = "e3w_cellcentered/mesh_mask.nc"
# """ e3w_watercellcentered """
# pdtB = "e3w_watercellcentered/overflow_zps_e3w_watercentered_grid_T.nc"
# pmmB = "e3w_watercellcentered/mesh_mask.nc"
# """ Diff """
# ptitle = r"Cellcentered - Watercell centered "
# mtitle = "diffT.mp4"


# """ e3w_cellcentered """
# pdtA = "e3w_cellcentered_e3u/overflow_zps_e3w_cellcentered_e3u_grid_T.nc"
# pmm = "e3w_cellcentered_e3u/mesh_mask.nc"
# """ e3w_watercellcentered """
# pdtB = "e3w_watercellcentered/overflow_zps_e3w_watercentered_grid_T.nc"
# pmmB = "e3w_watercellcentered/mesh_mask.nc"
# """ Diff """
# ptitle = r"Cellcentered_e3u - Watercell centered "
# mtitle = "diffT.mp4"

dtA = nc4.Dataset(pdtA)
thetaoA = dtA.variables['thetao_inst']
thetaoA = np.ma.masked_where(thetaoA[:]<5.,thetaoA[:])

dtB = nc4.Dataset(pdtB)
thetaoB = dtB.variables['thetao_inst']
thetaoB = np.ma.masked_where(thetaoB[:]<5.,thetaoB[:])

theta = thetaoA - thetaoB

""" Cell centered grid """

mm = nc4.Dataset(pmm)
nT,_,_,_ = np.shape(theta)
e3w = mm.variables['e3w_0'][0,:,1,:]
e1u = mm.variables['e1u'][0,1,:]

gridk = np.zeros((e3w.shape[0]+1, e3w.shape[1]+1))
gridk[1:,:-1] = np.cumsum(e3w,axis=0)
gridk[1:,-1]  = gridk[1:,-2]

gridx = np.zeros((e3w.shape[0]+1, e3w.shape[1]+1))
gridx[:,1:] = np.repeat([np.cumsum(e1u)], repeats = gridk.shape[0],axis = 0)

nK,_ = np.shape(gridk[1:,:-1])
gridp  = np.copy(gridk[1:,:-1])
gridxp = np.copy(gridx[1:,:-1])
for n in range(nK):
    if (n%3!=0):
        gridp[n,:] = np.nan
gridp = np.ma.masked_where(theta[0,:,1,:]<5.,gridp)

""" Movie
"""

vmin = -0.5 ; vmax = 0.5 ; dv = 0.1
levels = np.arange(vmin,vmax,dv)
levelc = np.arange(vmin,vmax,1)
palette = plt.get_cmap('bwr')

""" figure bordelum """
fig, ax = plt.subplots(figsize = [12, 8] ,dpi=200)

im = plt.pcolormesh(gridx, gridk, theta[0,:,1,:],
                    vmin=vmin, vmax=vmax,
                    cmap = palette)
cbar = plt.colorbar(im)
cbar.set_label(r"Temperature ($\degree C$)")
ax.plot(gridxp.T,gridp.T, 'k-', lw=0.5)

ax.set_ylim(2020,0)
ax.set_yticks([0,500,1000,1500,2000])
ax.set_yticklabels(["0","500","1000","1500","2000"])
ax.set_ylabel("depth (m)")

ax.set_xlim(0,200000)
ax.set_xticks([0,50000, 100000, 150000, 200000])
ax.set_xticklabels(["0","50","100","150","200"])
ax.set_xlabel("length (km)")

##################################################

ax.patch.set_color('1.')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", which = 'both', width=1.5, labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1.5, labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
plt.tight_layout(rect=[0,0,1,0.95])

""" movie thing """
def update_img(n,nT,ax,im,theta,ptitle):
    import sys
    sys.stdout.write(u"\u001b[1000D" + "processing [%3d/%3d]" % (n+1,nT))
    sys.stdout.flush()
    ax.set_title(ptitle + '%02d/%02d '%(n+1,nT), fontsize = 18, y = 1.02)
    im.set_array(theta[n,:,1,:].ravel())
    return

# https://stackoverflow.com/questions/18797175/animation-with-pcolormesh-routine-in-matplotlib-how-do-i-initialize-the-data
ani = animation.FuncAnimation(fig,update_img,  fargs = (nT,ax,im,theta,ptitle,) ,frames = nT, blit=False,repeat=False)
writer = animation.writers['ffmpeg'](fps=6)
ani.save(psave,writer=writer,dpi=200)
print("\nsaving : %s" % psave)
                                                                                                                                                                                                                                                                                                                      emAM.py                                                                                             000644  000766  000024  00000031074 14052143570 012024  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         # -*- coding: utf-8 -*-

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
#                     help="rotation angle of the grid (default 0??)")
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
# pdt = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_REF_16/freeslip_0/AM98_ref_16_freeslip_0_5d_00010101_00051230_grid_T.nc"
# pdu = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_REF_16/freeslip_0/AM98_ref_16_freeslip_0_10y_5d_00010101_00101230_grid_U.nc"
# pdv = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_REF_16/freeslip_0/AM98_ref_16_freeslip_0_10y_5d_00010101_00101230_grid_V.nc"
# pmm = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_REF_16/freeslip_0/mesh_mask.nc"
# dx = 6.25 ; dy =  6.25 # 25 km
# dx*=1000;dy*=1000
# theta = 0.
# vmin = -0.25 ; vmax = 0.25


pdt = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_ref4/freeslip0/AM98_ref_4_freeslip_0_5d_00010101_00051230_grid_T.nc"
pdu = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_ref4/freeslip0/AM98_ref_4_freeslip_0_5d_00010101_00051230_grid_U.nc"
pdv = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_ref4/freeslip0/AM98_ref_4_freeslip_0_5d_00010101_00051230_grid_V.nc"
pmm = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_ref4/freeslip0/mesh_mask.nc"
dx = 25 ; dy =  25 # 25 km
dx*=1000;dy*=1000
theta = 0.
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

# grid centr?? sur les corners bas gauche
lx = dx*np.cos((45+theta)*np.pi/180)/np.sqrt(2)
ly = dx*np.sin((45+theta)*np.pi/180)/np.sqrt(2)
# ballec cotes droit/haut - ils seront crop??s
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
# linear friction
wind = np.zeros((nI,nJ)) ; fric = np.zeros((nI,nJ)) ; lap = np.zeros((nI,nJ))
z2d = hu*uu*uu + hv*vv*vv
fric = - rho0 * rfr * z2d
# wind
z2d = tu*uu + tv*vv
wind = z2d
# laplacian
lap = rho0 * divgrad

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

# titlezer += "min = %3.1f m   max = %3.1f m   " % (np.min(data), np.max(data)) \
#          +  "$\Delta$ = %3.1f m" % (np.max(data) - np.min(data))
# ax.set_title(titlezer,
#               fontsize = 10, y = 1.02)
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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                    hello.py                                                                                            000755  000766  000024  00000000107 13762127137 012314  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         #!/usr/bin/env python3
import sys

for arg in sys.argv:
    print(arg)
                                                                                                                                                                                                                                                                                                                                                                                                                                                         plothistrhox.py                                                                                     000644  000766  000024  00000012307 14004626450 013752  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         #!/usr/bin/env python
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
# import cmocean
import xarray as xr


""" *****************************************************************
"""

import argparse
parser = argparse.ArgumentParser(description = "Plot the e3 weighted density histogram or overflow test case")
parser.add_argument("netfile", type=str,
                    help="Netcdf file path")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
parser.add_argument("-nt", type=float,
                    default=0.,
                    help="timeframe used")
parser.add_argument("-min","--minimum", type=float,
                    default = 10.,
                    help="minimal value of the colorbar")
parser.add_argument("-max","--maximum", type=float,
                    default = 20.,
                    help="maximal value of the colorbar")
parser.add_argument("-b", type=int,
                    default=50,
                    help="number of bins")
parser.add_argument("-s","--save", action="count", default=0,
                    help="save the current figure")
parser.add_argument("-m","--meshmask", type=str,
                    default="mesh_mask.nc",
                    help="meshmask associated to the file")

args = parser.parse_args()

if args.verbose:
    print("Running '{}'".format(__file__))
    print("netcdf used " + args.netfile)

""" *****************************************************************
"""
# nt = args.nt
pdt = args.netfile

# Result in, if there is no min/max declared, =None and pcolormesh handle this
Tmin = args.minimum ; Tmax = args.maximum

try:
    pmm = args.meshmask
    nc4.Dataset(pmm)
except:
    try:
        pmm = "/Users/gm/Documents/pythonGear/meshmaskOVERFLOW/"+args.meshmask
        nc4.Dataset(pmm)
    except:
        exit
""" *****************************************************************
"""

""" *****************************************************************
    Plot du test overflow - zps : z partial cells
    34 : time step
    101 : depth (z)
    3 : y dimension (1 cell embraced by two boundary layers)
    202 : x dimension
"""
# dt =  xr.open_dataset(pdt)
dt = nc4.Dataset(pdt)
mm = nc4.Dataset(pmm)

thetao = dt.variables['thetao_inst']
nT,_,_,nX = np.shape(thetao)
theta = np.ma.masked_outside(thetao, Tmin-0.5,Tmax)

e3t = dt.variables['e3t_inst']
e1u = mm.variables['e1u'][0,1,:]

gridx = np.cumsum(e1u)

""" figure bordelum """
fig, ax = plt.subplots(figsize=(8,4))

# https://stackoverflow.com/questions/21532667/numpy-histogram-cumulative-density-does-not-sum-to-1
# as the bins are 1/5 withd, the sum does equals to 1 (stepfilled)
kwargs = dict(histtype='bar', alpha=0.9, density=False, bins=args.b)


# # 30min time frame, 34 in total (17h)
# for t in [18,12,6,0]:
#     #        toce   x y z
#     data = (theta[t,:,1,:]    ).flatten()
#     e3   = (  e3t[t,:,1,:]/20.).flatten()
#     plt.hist(x=data,weights=e3,
#              label=r"%dh" % (t*0.5),
#              **kwargs)

# multiplebar
# 30min time frame, 34 in total (17h)
#                  # in km
Xframe = np.array([38,45,100])*1000 # in m
timeframe = []
for ix in range(len(Xframe)):
    il = np.argmin([np.abs(Xframe[ix] - gridx[i]) for i in range(nX)])
    t=0 ; flag = True
    while ( flag and t<(nT-1) ):
        test = np.any(theta[t,:,1,il]<20.)
        if test:
            flag=False
        else:
            t+=1
    # print("%3fkm" % (Xframe[il]/1000))
    # print("%3fkm" % (Xframe[il]/1000))
    timeframe.append(t)

data = [(theta[t,:,1,:]    ).flatten() for t in timeframe]
e3   = [(  e3t[t,:,1,:]/20.).flatten() for t in timeframe]
plt.hist(x=data,weights=e3,
         label=[r"%3.0fkm (%.1fh)" %   (Xframe[it[0]]/1E3, it[1]*0.5)
                                    for it in zip(range(len(Xframe)), timeframe)],
          **kwargs)

titlezer = "%s\n" % pdt
titlezer += r"e3t Weighted $\theta$ Histogram (total of %dh)" % (nT*0.5)
psave = "histrhox.png"

ax.set_title(titlezer,y = 1.02)
#
ax.set_xlim((Tmin,Tmax))
ax.set_xlabel(r"$\theta$ (??C)")

ax.set_ylim(0,600)
ax.set_ylabel("Density")
ax.grid(which='major', linestyle='-', linewidth='0.5', color='red')
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')

plt.legend(loc=9, ncol=3)
##################################################

ax.patch.set_color('1.')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", which = 'both', width=1.5, labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1.5, labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
plt.tight_layout(rect=[0,0,1,0.95])
plt.tight_layout()

# showsave(psave, fig, save)
if args.save :
    print("saving : %s" % psave)
    fig.savefig(psave)
else :
    plt.show()
                                                                                                                                                                                                                                                                                                                         qRidge.py                                                                                           000644  000766  000024  00000010526 14030346324 012415  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         #!/usr/bin/env python
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
import gsw


""" *****************************************************************
"""

import argparse
parser = argparse.ArgumentParser(description = "Plot the down flow of a denser water mass (zcoor) with the overflow config")
parser.add_argument("netfile", type=str,
                    help="Netcdf U grid file path")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
parser.add_argument("-s","--save", action="count", default=0,
                    help="save the current figure")
parser.add_argument("-t","--pdt", type=str,
                    default="",
                    help="netfile on T point")
parser.add_argument("-m","--meshmask", type=str,
                    default="mesh_mask.nc",
                    help="meshmask associated to the file")

args = parser.parse_args()

if args.verbose:
    print("Running '{}'".format(__file__))
    print("netcdf used " + args.netfile)
########################################################

pdu = args.netfile
pdt = args.pdt

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
psave = "qridge.png"
########################################################


# pdt = "ridge/sco/RIDGE_ref_sco_12h_grid_T.nc"
# pdu = "ridge/sco/RIDGE_ref_sco_12h_grid_U.nc"
# pmm = "ridge/sco/mesh_mask.nc"
# save = 0

dtu = nc4.Dataset(pdu)
dt  = nc4.Dataset(pdt)

uu  = dtu.variables['u_vol']
toce = dt.variables['toce']
nT,nK,nY,nX = np.shape(uu)

mm = nc4.Dataset(pmm)
umask = mm.variables['umask']
tmask = mm.variables['tmask']

toce = np.ma.masked_where(tmask==0,toce)
uum  = np.ma.masked_where(umask==0,uu  )

glam = dt.variables['nav_lon'] ; gphi = dt.variables['nav_lat']

midY = nY//2+1
midX = nX//2+1

def tfromrho(rho):
    f = gsw.CT_from_rho(1000+rho, 34.7, 4000)[0]
    f = gsw.t_from_CT(34.7, f, 4000)
    return(f)

Tlim = np.array([45.,45.5,46.,46.5])
nR = np.shape(Tlim)[0]
Q = np.zeros((nT, nR-1))

timelist = (np.arange(0,nT,1)+1.)/2.
""" transport per classe """
for t in range(nT):
    for r in range(nR-1):
        # r = 0
        ta=Tlim[r] ; tb = [r+1]
        # classe de temperature
        tmpm= np.ma.masked_outside(toce[t,:,:,midX],ta,tb)
        # flow sortant de la classe de temperature
        tmp = np.ma.masked_where(tmpm.mask, uum[t,:,:,midX])
        tmp = np.ma.masked_where(tmp<0, tmp) # la somme globale abs<0.01 Sv
        q = np.nansum(np.ma.filled(tmp,0.))
        Q[t,r] = q/1E6

""" figure """
fig, ax = plt.subplots(figsize=(8,4), dpi=200)
# at each front
for r in range(nR-1):
    ax.plot(timelist, Q[:,r],
            label=r"$%2.2f<\rho<%2.2f$ [kg/m3] | $%2.2f>T>%2.2f$ [??C]" % \
            ( Tlim[r],Tlim[r+1],tfromrho(Tlim[r]),tfromrho(Tlim[r+1]) ),
            # label=r"$%2.2f<\rho<%2.2f$ [kg/m3]" % \
            # ( Tlim[r], Tlim[r+1] ),
             alpha=1., linewidth=.7)
    # conversion densit??/temp??rature


ax.grid(which='major', linestyle='-', linewidth='0.3', color='black')
ax.grid(which='minor', linestyle=':', linewidth='0.3', color='black')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
ax.set_xlim(timelist[0],timelist[-1])
ax.set_xlabel("Time in days")
ax.set_ylim(0,8)
ax.set_ylabel("Volume Transport (Sv)")
ax.legend(loc=9, ncol=1)

fig.add_subplot(ax)

# print(toce[-1,:,4,4])
titlezer  = '%s\n'%(pdt)
titlezer += "min = %2.2f kg/m3   max = %2.2f kg/m3   " % (np.min(toce[toce!=0]), np.max(toce[toce!=0]))
ax.set_title(titlezer, fontsize = 12, y = 1.02)
plt.tight_layout(rect=[0,0,1,0.95])


if save :
    print("saving : %s" % psave)
    fig.savefig(psave)
else :
    plt.show()
                                                                                                                                                                          rhoOverflow.py                                                                                      000644  000766  000024  00000022631 14051425055 013520  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         #!/usr/bin/env python
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
parser = argparse.ArgumentParser(description = "Plot the e3 weighted density histogram or overflow test case")
parser.add_argument('-n','--netfile', nargs='+',
                    type=str, required=True,
                    help="Netcdf file path")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
parser.add_argument("-eps","--epsilon", type=float,
                    default = 1E-2,
                    help="tolerance on temperature")
parser.add_argument("-dT", type=float,
                    default=0.5,
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
Tbins = np.arange(10.,20.,dT)
# Tbins = np.append(Tbins,20.)
""" Front position """
# Xframe = np.array([40,60,100])*1000 # in m
# Xframe = np.array([60,100])*1000 # in m
Xframe = np.array([60])*1000 # in m

colors=["royalblue","orangered","forestgreen", "black"]
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

    """ corresponding time frame """
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
    to bring them back to 10?? volume
    """
    Vcold = 0
    for the in range(len(Tbins)):
        # print("[%1.2f,%1.2f[" % ((Tbins[the]-dT/2.)-eps, (Tbins[the]+dT/2.)) )
        tmp = np.ma.masked_outside(thetao[0,:,1,:],
                                   (Tbins[the]-dT/2.)-eps, (Tbins[the]+dT/2.))
        # volume d'eau ?? the
        e3 = np.ma.masked_where(tmp.mask, e3t[0,:,1,:])
        e3 = np.nansum(np.ma.filled(e3,0.))
        tmp2 = e3
        tmp2 *= (20.-Tbins[the])/10.
        Vcold+=tmp2
    # print("%f" % Vcold)
    """ Distribution """
    for t in range(len(Xframe)):
        for the in range(len(Tbins)):
            tmp = np.ma.masked_outside(thetao[timeframe[ndt,t],:,1,:],
                                       (Tbins[the]-dT/2.)-eps, (Tbins[the]+dT/2.))
            # volume d'eau ?? the
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
# ax[1].set_xlabel("Temperature (??C)")
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
                marker=markers[t], color=colors[ndt],
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
ax.set_xlabel("Temperature (??C)")
ax.set_ylim(0,0.3)
ax.set_ylabel("Dense water corresponding volume ratio")
ax.legend(loc=9, ncol=len(Xframe))

fig.add_subplot(ax)

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
                                                                                                       rhoRidge.py                                                                                         000644  000766  000024  00000022662 14026373226 012757  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         #!/usr/bin/env python
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
    to bring them back to 10?? volume
    """
    Vcold = 0
    for the in range(len(Tbins)):
        # print("[%1.2f,%1.2f[" % ((Tbins[the]-dT/2.)-eps, (Tbins[the]+dT/2.)) )
        tmp = np.ma.masked_outside(thetao[0,:,1,:],
                                   (Tbins[the]-dT/2.)-eps, (Tbins[the]+dT/2.))
        # volume d'eau ?? the
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
            # volume d'eau ?? the
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
# ax[1].set_xlabel("Temperature (??C)")
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
ax.set_xlabel("Temperature (??C)")
ax.set_ylim(0,0.3)
ax.set_ylabel("Dense water corresponding volume ratio")
ax.legend(loc=9, ncol=len(Xframe))

fig.add_subplot(ax)

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
                                                                              sfRidge.py                                                                                          000644  000766  000024  00000021064 14030352762 012567  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         # -*- coding: utf-8 -*-

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

pdt = "ridge/sco/RIDGE_ref_sco_12h_grid_T.nc"
pmm = "ridge/sco/mesh_mask.nc"
pdu = "ridge/sco/RIDGE_ref_sco_12h_grid_U.nc"

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

uu  = dtu.variables['u_vol']   # volumic transport
toce = dt.variables['toce']
nT,nK,nY,nX = np.shape(uu)

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

# gridk et gridx sont les sommets d??finissant les pixels de pcolormesh
gridk = np.zeros((e3t.shape[0]+1, e3t.shape[1]+1))
gridk[1:,:-1] = np.cumsum(e3t,axis=0) # somme sur la vertical de i,j=0,0 (upper left corner)
                                      # vers les i ascendant (z descendant)
gridk[1:,-1]  = gridk[1:,-2]          # duplique la derni??re colonne avec l'avant derni??re

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
            sfm[t,k,i] += - np.sum(uum[t,k:,:,i])/1E6

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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                            sshAM.py                                                                                            000644  000766  000024  00000016162 14051476276 012234  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         # -*- coding: utf-8 -*-

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

import argparse
parser = argparse.ArgumentParser(description = "Plot the ssh for AM98")
parser.add_argument("netfile", type=str,
                    help="Netcdf file path")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")

parser.add_argument("-dx", type=real,
                    default=25,
                    help="grid spacing in km")
parser.add_argument("-t","--theta", type=real,
                    default=0.,
                    help="rotation angle of the grid (default 0??)")
parser.add_argument("-s","--save", action="count", default=0,
                    help="save the current figure")
parser.add_argument("-f","--film", action="count", default=0,
                    help="make a film")

parser.add_argument("-min","--minimum", type=float, default = 200,
                    help="minimal value of the colorbar")
parser.add_argument("-max","--maximum", type=float, default = 800,
                    help="maximal value of the colorbar")

parser.add_argument("-w","--width", type=float, default = 0.,
                    help="should the band be reprenseted in the graph")

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
dx =  float(args.dx) ; dy =  float(args.dx) # 25 km
dx*=1000;dy*=1000
pdtT = args.netfile

print(args.meshmask)
# Result in, if there is no min/max declared, =None and pcolormesh handle this
vmin = args.minimum ; vmax = args.maximum

try:
    pmm = args.meshmask
    nc4.Dataset(pmm)
except:
    try:
        pmm = "../"+args.meshmask
        nc4.Dataset(pmm)
    except:
        try:
            pmm = "/Users/gm/Documents/pythonGear/meshmaskAM98/"+args.meshmask
            nc4.Dataset(pmm)
        except:
            exit
""" *****************************************************************
"""

# ssh
dt = nc4.Dataset(pdtT)

try:
    ssh = dt.variables['sossheig'][:,:,:]
except:
    ssh = dt.variables['ht'][:,:,:]
nT,_,_ = np.shape(ssh)
mm = nc4.Dataset(pmm)
tmask = mm.variables['tmask'][0,0]

# glamt, gphit
try:
    glamt = dt.variables['nav_lon'][:]
    gphit = dt.variables['nav_lat'][:]
except:
    glamt = mm.variables['nav_lon'][:]
    gphit = mm.variables['nav_lat'][:]

# grid centr?? sur les corners bas gauche
lx = dx*np.cos((45+theta)*np.pi/180)/np.sqrt(2)
ly = dx*np.sin((45+theta)*np.pi/180)/np.sqrt(2)
# ballec cotes droit/haut - ils seront crop??s
nx = np.shape(glamt)[0] ; ny = np.shape(glamt)[1]
gridx = np.zeros((nx,ny))
gridx = glamt - lx
gridy = np.zeros((nx,ny))
gridy = gphit - ly

limx = dx ; limy = dx
# limx = dx ; limy = dx
# limx = dx*np.cos((45+theta)*np.pi/180)
# limy = dx*np.sin((45+theta)*np.pi/180)

ssh=ssh+500

""" SSH map (end experiment)
"""

########################################################
if args.film :
    data = ssh[0,:,:]
    data = np.ma.masked_where(tmask==0,data)
    titlezer  = '%s - %02d/%02d \n'%(pdtT,0+1,nT)
    psave = "ssh.mp4"
else :
    data = ssh[-1,:,:]
    data = np.ma.masked_where(tmask==0,data)
    titlezer  = '%s\n'%(pdtT)
    psave = "ssh.png"
########################################################

palette = plt.get_cmap('RdYlBu_r')
# palette = plt.get_cmap('seismic')

""" figure bordelum """
fig, ax = plt.subplots(dpi=200)

im = ax.pcolormesh(gridx, gridy, data,
                    vmin=vmin, vmax=vmax,
                    # levels = levels,
                    cmap = palette)
cbar = plt.colorbar(im)
if args.film==0 :
    c1 = ax.contour(glamt,gphit, data,
                    # levels = levelc,
                    vmin=vmin,vmax=vmax,
                    linewidths =0.3, colors=('k',),linestyles = "solid")
    ax.clabel(c1, fmt='%3.0f', colors='k', fontsize=10)

cbar.set_label(r"Upper Layer Width - (m)")

if args.width == 0.:
    ax.set_xlim(-limx,2000000)
    ax.set_xticks([0,500000,1000000,1500000,2000000])
    ax.set_xticklabels(["0","500","1000","1500","2000"])
    ax.set_xlabel("X (km)")

    ax.set_ylim(-limy,2000000)
    ax.set_yticks([0,500000,1000000,1500000,2000000])
    ax.set_yticklabels(["0","500","1000","1500","2000"])
    ax.set_ylabel("Y (km)")

titlezer += "min = %3.1f m   max = %3.1f m   " % (np.min(data), np.max(data)) \
         +  "$\Delta$ = %3.1f m" % (np.max(data) - np.min(data))
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


if args.film :
    """ movie thing """
    def update_img(n,nT,ax,im,theta,ptitle):
        import sys
        data = ssh[n,:,:]
        data = np.ma.masked_where(tmask==0,data)
        ptitle = '%s - %02d/%02d \n'%(pdtT,n+1,nT)
        ptitle += "min = %3.1f m   max = %3.1f m   " % (np.min(data), np.max(data)) \
               +  "$\Delta$ = %3.1f m" % (np.max(data) - np.min(data))
        sys.stdout.write(u"\u001b[1000D" + "processing [%3d/%3d]" % (n+1,nT))
        sys.stdout.flush()
        ax.set_title(ptitle, fontsize = 10, y = 1.02)
        # new_color = palette(data.T.ravel())
        #
        # for tp in c1.collections:
        #     tp.remove()
        # c1 = ax.contour(glamt,gphit, data,
        #                 # levels = levelc,
        #                 vmin=vmin,vmax=vmax,
        #                 linewidths =0.3, colors=('k',),linestyles = "solid")
        # ax.clabel(c1, fmt='%3.0f', colors='k', fontsize=10)
        #
        # im.update({'facecolors':new_color})
        data = data[:-1, :-1]
        # ravel() converts C to a 1d-array
        im.set_array(data.ravel())
        return

    # https://stackoverflow.com/questions/18797175/animation-with-pcolormesh-routine-in-matplotlib-how-do-i-initialize-the-data
    ani = animation.FuncAnimation(fig,update_img,  fargs = (nT,ax,im,ssh,titlezer,) ,frames = nT, blit=False,repeat=False)
    writer = animation.writers['ffmpeg'](fps=nT/10)
    ani.save(psave,writer=writer,dpi=200)
    print("\nsaving : %s" % psave)
else :
    if args.save :
        print("\nsaving : %s" % psave)
        fig.savefig(psave)
        plt.close()
    plt.show()
                                                                                                                                                                                                                                                                                                                                                                                                              test.py                                                                                             000644  000766  000024  00000005052 13760752170 012170  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         # -*- coding: utf-8 -*-

# import argparse
# # https://docs.python.org/3/library/argparse.html
# parser = argparse.ArgumentParser(description='Add some integers.')
#
# # argument positionnel
# parser.add_argument('integers', metavar='N', type=int,
#                     nargs='+', # iterable
#                     help='integer list')
# # argument optionnel
# # action='store_const' = si rien n'est pr??cis??, active tout de m??me ?
# # dest : nom qu'on lui donne, autrement il s'appelerait sum
# parser.add_argument('-s','--sum', dest='accumulate', action='store_const',
#                     const=sum, default=max,
#                     help='sum the integers (default: find the max)')
#
#
# args = parser.parse_args()
# print(args.accumulate(args.integers))

""" https://docs.python.org/fr/3/howto/argparse.html
    Premier exemple parlant
"""
# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument("square", type=int,
#                     help="display a square of a given number")
# parser.add_argument("-v", "--verbose", action="store_true",
#                     help="increase output verbosity")
# args = parser.parse_args()
# answer = args.square**2
# if args.verbose:
#     print("the square of {} equals {}".format(args.square, answer))
# else:
#     print(answer)

""" https://docs.python.org/fr/3/howto/argparse.html
    Deuxi??me exemple, plus pratique
"""
# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument("x", type=int, help="the base")
# parser.add_argument("y", type=int, help="the exponent")
# parser.add_argument("-v", "--verbosity", action="count", default=0)
# args = parser.parse_args()
# answer = args.x**args.y
# if args.verbosity >= 2:
#     print("Running '{}'".format(__file__))
# if args.verbosity >= 1:
#     print("{}^{} == ".format(args.x, args.y), end="")
# print(answer)

""" https://docs.python.org/fr/3/howto/argparse.html
    Troisi??me exemple conflit
"""
# import argparse
#
# parser = argparse.ArgumentParser(description="calculate X to the power of Y")
# group = parser.add_mutually_exclusive_group()
# group.add_argument("-v", "--verbose", action="store_true")
# group.add_argument("-q", "--quiet", action="store_true")
# parser.add_argument("x", type=int, help="the base")
# parser.add_argument("y", type=int, help="the exponent")
# args = parser.parse_args()
# answer = args.x**args.y
#
# if args.quiet:
#     print(answer)
# elif args.verbose:
#     print("{} to the power {} equals {}".format(args.x, args.y, answer))
# else:
#     print("{}^{} == {}".format(args.x, args.y, answer))

""" Mon exemple
"""
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      toceOverflow.py                                                                                     000755  000766  000024  00000014320 14027137136 013664  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         # -*- coding: utf-8 -*-

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
parser = argparse.ArgumentParser(description = "Plot the down flow of a denser water mass (zcoor) with the overflow config")
parser.add_argument("netfile", type=str,
                    help="Netcdf file path")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
parser.add_argument("-s","--save", action="count", default=0,
                    help="save the current figure")

parser.add_argument("-nt", type=int,
                    default=-1,
                    help="timeframe used")
parser.add_argument("-f","--film", action="count", default=0,
                    help="make a film")

parser.add_argument("-dv", type=float, default = 0.5,
                    help="color interval in the colorbar")
parser.add_argument("-min","--minimum", type=float, default = 10.,
                    help="minimal value of the colorbar")
parser.add_argument("-max","--maximum", type=float, default = 20.,
                    help="maximal value of the colorbar")

parser.add_argument("-m","--meshmask", type=str,
                    default="mesh_mask.nc",
                    help="meshmask associated to the file")

args = parser.parse_args()

if args.verbose:
    print("Running '{}'".format(__file__))
    print("netcdf used " + args.netfile)
########################################################

nt = args.nt
pdt = args.netfile
# Result in, if there is no min/max declared, =None and pcolormesh handle this
vmin = args.minimum ; vmax = args.maximum
try:
    pmm = args.meshmask
    nc4.Dataset(pmm)
except:
    try:
        pmm = "/Users/gm/Documents/pythonGear/meshmaskOVERFLOW/"+args.meshmask
        nc4.Dataset(pmm)
    except:
        exit

if args.dv==-1.:
    dv = (vmax-vmin)/100.
else:
    dv = args.dv

Ncolor=int((vmax-vmin)/dv)

########################################################

"""
    Plot du test overflow - zps : z partial cells
    34 : time step
    101 : depth (z)
    3 : y dimension (1 cell embraced by two boundary layers)
    202 : x dimension
"""
dt = nc4.Dataset(pdt)
mm = nc4.Dataset(pmm)

thetao = dt.variables['thetao_inst']
nT,_,_,_ = np.shape(thetao)
theta = np.ma.masked_where(thetao[:]<5.,thetao[:])

e3w = mm.variables['e3w_0'][0,:,1,:]
e1u = mm.variables['e1u'][0,1,:]

# gridk et gridx sont les sommets d??finissant les pixels de pcolormesh
gridk = np.zeros((e3w.shape[0]+1, e3w.shape[1]+1))
gridk[1:,:-1] = np.cumsum(e3w,axis=0) # somme sur la vertical de i,j=0,0 (upper left corner)
                                      # vers les i ascendant (z descendant)
gridk[1:,-1]  = gridk[1:,-2]          # duplique la derni??re colonne avec l'avant derni??re

gridx = np.zeros((e3w.shape[0]+1, e3w.shape[1]+1))
gridx[:,1:] = np.repeat([np.cumsum(e1u)], repeats = gridk.shape[0],axis = 0)

nK,_ = np.shape(gridk[1:,:-1])
gridp  = np.copy(gridk[1:,:-1])
gridxp = np.copy(gridx[1:,:-1])
for n in range(nK):
    if (n%3!=0):
        gridp[n,:] = np.nan
gridp = np.ma.masked_where(thetao[0,:,1,:]<5.,gridp)

########################################################
if args.film :
    data = theta[0,:,1,:]
    titlezer  = '%s - %02d/%02d \n'%(pdt,0+1,nT)
    psave = "overflow.mp4"
else :
    data = theta[nt,:,1,:]
    titlezer  = '%s - %02d/%02d \n'%(pdt,nt+1 if nt!=-1 else nT,nT)
    psave = "overflow.png"
########################################################
# levelc = np.arange(vmin,vmax,1)
palette = plt.get_cmap('jet',Ncolor)

""" figure bordelum """
fig, ax = plt.subplots(figsize = [12, 8])

im = plt.pcolormesh(gridx, gridk, data,
                    cmap = palette, vmin=vmin, vmax=vmax)
cbar = plt.colorbar(im)
cbar.set_label(r"Temperature ($\degree C$)")
ax.plot(gridxp.T,gridp.T, 'k-', lw=0.5)

titlezer += "min = %3.1f ??C   max = %3.1f ??C   " % (np.min(data), np.max(data)) +  "$\Delta$ = %3.1f ??C" % (np.max(data) - np.min(data))
ax.set_title(titlezer, fontsize = 18, y = 1.02)

ax.set_ylim(2020,0)
ax.set_yticks([0,500,1000,1500,2000])
ax.set_yticklabels(["0","500","1000","1500","2000"])
ax.set_ylabel("depth (m)")

ax.set_xlim(0,200000)
ax.set_xticks([0,50000, 100000, 150000, 200000])
ax.set_xticklabels(["0","50","100","150","200"])
ax.set_xlabel("length (km)")

ax.patch.set_color('1.')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
plt.tight_layout(rect=[0,0,1,0.95])
# plt.tight_layout()

##################################################

if args.film :
    """ movie thing """
    def update_img(n,nT,ax,im,theta,ptitle):
        import sys
        data = theta[n,:,1,:]
        ptitle = '%s - %02d/%02d \n'%(pdt,n+1,nT)
        ptitle += "min = %3.1f ??C   max = %3.1f ??C   " % (np.min(data), np.max(data)) +  "$\Delta$ = %3.1f ??C" % (np.max(data) - np.min(data))

        sys.stdout.write(u"\u001b[1000D" + "processing [%3d/%3d]" % (n+1,nT))
        sys.stdout.flush()
        ax.set_title(ptitle, fontsize = 18, y = 1.02)
        im.set_array(data.ravel())
        return

    # https://stackoverflow.com/questions/18797175/animation-with-pcolormesh-routine-in-matplotlib-how-do-i-initialize-the-data
    ani = animation.FuncAnimation(fig,update_img,  fargs = (nT,ax,im,theta,titlezer,) ,frames = nT, blit=False,repeat=False)
    writer = animation.writers['ffmpeg'](fps=6)
    ani.save(psave,writer=writer,dpi=200)
    print("\nsaving : %s" % psave)
else :
    if args.save :
        print("\nsaving : %s" % psave)
        fig.savefig(psave)
    plt.show()
##################################################
                                                                                                                                                                                                                                                                                                                toceRidge.py                                                                                        000755  000766  000024  00000017156 14027170672 013127  0                                                                                                    ustar 00gm                              staff                           000000  000000                                                                                                                                                                         # -*- coding: utf-8 -*-

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
parser = argparse.ArgumentParser(description = "Plot the down flow of a denser water mass (zcoor) with the rige config")
parser.add_argument("netfile", type=str,
                    help="Netcdf file path")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
parser.add_argument("-s","--save", action="count", default=0,
                    help="save the current figure")

parser.add_argument("-nt", type=int,
                    default=-1,
                    help="timeframe used")
parser.add_argument("-f","--film", action="count", default=0,
                    help="make a film")

parser.add_argument("-dv", type=float, default = 0.02,
                    help="color interval in the colorbar")
parser.add_argument("-min","--minimum", type=float, default = 45.7,
                    help="minimal value of the colorbar")
parser.add_argument("-max","--maximum", type=float, default = 46.3,
                    help="maximal value of the colorbar")

parser.add_argument("-m","--meshmask", type=str,
                    default="mesh_mask.nc",
                    help="meshmask associated to the file")

args = parser.parse_args()

if args.verbose:
    print("Running '{}'".format(__file__))
    print("netcdf used " + args.netfile)
########################################################

nt = args.nt
pdt = args.netfile
# Result in, if there is no min/max declared, =None and pcolormesh handle this
vmin = args.minimum ; vmax = args.maximum
try:
    pmm = args.meshmask
    nc4.Dataset(pmm)
except:
    try:
        pmm = "/Users/gm/Documents/pythonGear/meshmaskOVERFLOW/"+args.meshmask
        nc4.Dataset(pmm)
    except:
        exit

if args.dv==-1.:
    dv = (vmax-vmin)/100.
else:
    dv = args.dv



########################################################


# pdt = "ridge/RIDGE_ref_zco_12h_grid_T.nc"
# pmm = "ridge/mesh_mask.nc"
# vmin = 45.7 ; vmax = 46.3 ; dv = 0.02

Ncolor=int((vmax-vmin)/dv)

"""
    Plot du test overflow - zps : z partial cells
    34 : time step
    101 : depth (z)
    3 : y dimension (1 cell embraced by two boundary layers)
    202 : x dimension
"""
dt = nc4.Dataset(pdt)
mm = nc4.Dataset(pmm)

thetao = dt.variables['toce']
nT,nK,nY,nX = np.shape(thetao)

print("domain size is (x,y) %dx%d with %d k-levels" % (nX,nY,nK))

middle = int(np.shape(thetao)[-2] / 2)
print("middle is %d" % (middle))
# middle = 3

tmask = mm.variables['tmask'][0,:,:,:]

sill=tmask[:,middle,:]
gauss=tmask[:,1,:]

e3w = mm.variables['e3w_0'][0,:,middle,:]
e3t = mm.variables['e3w_0'][0,:,middle,:]
                   #    time, x, y
# e1u = mm.variables['e1u'][0,middle,:]
# linx = dt.variables['nav_lon'][0,:]

e1t   = mm.variables['e1t'  ][0,middle,:]
glamt = mm.variables['glamt'][0,middle,:]

########################################################

# gridk et gridx sont les sommets d??finissant les pixels de pcolormesh
gridk = np.zeros((e3t.shape[0]+1, e3t.shape[1]+1))
gridk[1:,:-1] = np.cumsum(e3t,axis=0) # somme sur la vertical de i,j=0,0 (upper left corner)
                                      # vers les i ascendant (z descendant)
gridk[1:,-1]  = gridk[1:,-2]          # duplique la derni??re colonne avec l'avant derni??re

gridx = np.zeros(gridk.shape)       # grille x
# gridx[:,:-1] = linx - e1u/2         # sur toute la colonne,
# gridx[:, -1] = linx[-1]+e1u[-1]/2

gridx[:,:-1] = glamt    - e1t    /2         # sur toute la colonne,
gridx[:, -1] = glamt[-1]+ e1t[-1]/2

# nK,_ = np.shape(gridk)
gridkp = np.copy(gridk)
gridxp = np.copy(gridx)
for n in range(nK):
    if (n%2!=0):
        gridkp[n,:] = np.nan
# gridkp = np.ma.masked_where(thetao[0,:,middle,:]<1.,gridkp)

########################################################

# theta = np.ma.masked_where(thetao[:]<1.,thetao)
theta = np.copy(thetao)

ttmask = np.copy(thetao)
for t in range(nT):
    ttmask[t] = tmask
theta = np.ma.masked_where(ttmask<1,thetao)

for t in range(nT):
    theta[t,:,:,:] = np.ma.masked_where(tmask<1,thetao[t,:,:,:])
########################################################
if args.film :
    data = theta[0,:,middle,:]
    titlezer  = '%s - %02d/%02d \n'%(pdt,0+1,nT)
    psave = "ridge.mp4"
else :
    data = theta[nt,:,middle,:]
    titlezer  = '%s - %02d/%02d \n'%(pdt,nt+1 if nt!=-1 else nT,nT)
    psave = "ridge.png"

# titlezer  = '%s\n'%(pdt)
# data = theta[0,:,middle,:]

########################################################
# levelc = np.arange(vmin,vmax,1)
palette = plt.get_cmap('jet_r',Ncolor)

""" figure bordelum """
fig, ax = plt.subplots(figsize = [12, 8])

# plt.pcolormesh(gridx, gridk, gauss,alpha=1.,zorder = 0.)
im = plt.pcolormesh(gridx, gridk, data,
                    cmap = palette, vmin=vmin, vmax=vmax)
cbar = plt.colorbar(im)
cbar.set_label(r"Density Anomaly $\sigma_4$ (kg/m3)")
ax.plot(gridxp.T,gridkp.T, 'k-', lw=0.5)
# https://stackoverflow.com/questions/31877353/overlay-an-image-segmentation-with-numpy-and-matplotlib
# plt.pcolormesh(gridx, gridk, sill ,alpha=0.9,zorder = 0.)
# plt.pcolormesh(gridx, gridk, gauss,alpha=0.5,zorder = 0.)

titlezer += "min = %2.2f kg/m3   max = %2.2f kg/m3   " % (np.min(data), np.max(data))
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

ax.patch.set_color('1.')
# ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
plt.tight_layout(rect=[0,0,1,0.95])
# plt.tight_layout()

##################################################

if args.film :
    """ movie thing """
    def update_img(n,nT,ax,im,theta,ptitle):
        import sys
        data = theta[n,:,middle,:]
        ptitle = '%s - %02d/%02d \n'%(pdt,n+1,nT)
        ptitle += "min = %2.2f ??C   max = %2.2f ??C   " % (np.min(data), np.max(data)) +  "$\Delta$ = %2.2f ??C" % (np.max(data) - np.min(data))

        sys.stdout.write(u"\u001b[1000D" + "processing [%3d/%3d]" % (n+1,nT))
        sys.stdout.flush()
        ax.set_title(ptitle, fontsize = 18, y = 1.02)
        im.set_array(data.ravel())
        return

    # https://stackoverflow.com/questions/18797175/animation-with-pcolormesh-routine-in-matplotlib-how-do-i-initialize-the-data
    ani = animation.FuncAnimation(fig,update_img,  fargs = (nT,ax,im,theta,titlezer,) ,frames = nT, blit=False,repeat=False)
    writer = animation.writers['ffmpeg'](fps=6)
    ani.save(psave,writer=writer,dpi=200)
    print("\nsaving : %s" % psave)
else :
    if args.save :
        print("\nsaving : %s" % psave)
        fig.savefig(psave)
    plt.show()
##################################################
# plt.show()
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  