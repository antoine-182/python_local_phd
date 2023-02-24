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
import gsw


""" *****************************************************************
"""

# import argparse
# parser = argparse.ArgumentParser(description = "Plot the transport (Sv) accros the ridge")
# parser.add_argument("netfile", type=str,
#                     help="Netcdf U grid file path")
# parser.add_argument("-v", "--verbose", action="store_true",
#                     help="increase output verbosity")
# parser.add_argument("-s","--save", action="count", default=0,
#                     help="save the current figure")
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
# pdu = args.netfile
# pdt = args.pdt
#
# try:
#     pmm = args.meshmask
#     nc4.Dataset(pmm)
# except:
#     try:
#         pmm = "/Users/gm/Documents/pythonGear/meshmaskOVERFLOW/"+args.meshmask
#         nc4.Dataset(pmm)
#     except:
#         exit

# save = args.save
psave = "qridge.png"
########################################################

dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref/"
pdt = dir+"gen_sco_dx10/RIDGE_12h_grid_T.nc"
pdu = dir+"gen_sco_dx10/RIDGE_12h_grid_U.nc"
pmm = dir+"gen_sco_dx10/mesh_mask.nc"
save = 0

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

midY = nY//2
midX = nX//2

def tfromrho(rho):
    f = gsw.CT_from_rho(1000+rho, 34.7, 4000)[0]
    f = gsw.t_from_CT(34.7, f, 4000)
    return(f)

# Tlim = np.array([45.,45.5,46.,46.5])
Tlim = np.array([45.,46.5])
nR = np.shape(Tlim)[0]
Q = np.zeros((nT, nR-1))

timelist = (np.arange(0,nT,1)+1.)/2. # timeframe every 12h
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
            # label=r"$%2.2f<\sigma_4<%2.2f$ [kg/m3] | $%2.2f>T>%2.2f$ [°C]" % \
            # ( Tlim[r],Tlim[r+1],tfromrho(Tlim[r]),tfromrho(Tlim[r+1]) ),
            # label=r"$%2.2f<\sigma_4<%2.2f$ [kg/m3]" % \
            # ( Tlim[r],Tlim[r+1] ),
            # label=r"$%2.2f<\rho<%2.2f$ [kg/m3]" % \
            # ( Tlim[r], Tlim[r+1] ),
            label ="_no_legend_",
             alpha=1., linewidth=.7)
    # conversion densité/température


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
# ax.set_ylim(0,8)
ax.set_ylabel("Volume Transport (Sv)")
# ax.legend(loc=9, ncol=1)

fig.add_subplot(ax)

# print(toce[-1,:,4,4])
# titlezer  = '%s\n'%(pdt)
titlezer = "Transport across the ridge\n"
# titlezer += "min = %2.2f kg/m3   max = %2.2f kg/m3   " % (np.min(toce[toce!=0]), np.max(toce[toce!=0]))
ax.set_title(titlezer, fontsize = 12, y = 1.02)
plt.tight_layout(rect=[0,0,1,0.95])


if save :
    print("saving : %s" % psave)
    fig.savefig(psave)
else :
    plt.show()
