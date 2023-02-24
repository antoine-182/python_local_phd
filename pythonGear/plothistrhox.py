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
ax.set_xlabel(r"$\theta$ (°C)")

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
