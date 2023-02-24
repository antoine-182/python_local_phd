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
