
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

dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_ref/ref3ubs"
dir = "/Users/gm/Documents/nemo/dev_14237_KERNEL-01_IMMERSE_SEAMOUNT/tests/ridge/EXP_bvp"
pmm = "/mesh_mask.nc"

rn_h_ridge = 2800. ;  rn_h_strait = 4500. ; S0 = 30*(rn_h_strait-rn_h_ridge)/1E3

mm = nc4.Dataset(dir+pmm)
tmask = mm.variables['tmask'][0] ; umask = mm.variables['umask'][0]
nK,nY,nX = np.shape(tmask)
glamt = mm.variables['glamt'][0] ; gphiv = mm.variables['gphiv'][0]
midY = np.where(np.abs(gphiv[:,0])<=1E3)[0][0]
midX = np.where(np.abs(glamt[0,:])<=1E3)[0][0]

e2u = mm.variables['e2u'][0,:,midX] ; e3u_0 = mm.variables['e3u_0'][0,:,:,midX] # penalise
gdepw = mm.variables['gdepw_0'][0]

a = umask[:,:,midX] ; b = gdepw[:,:,midX]
zone = (a==1)*(b>rn_h_ridge)  # sous le ridge
#
e2u_3D = np.repeat(e2u[np.newaxis,:],nK,axis=0)
S = np.sum(e2u_3D[zone]*e3u_0[zone])/1E6 # km2

print("Surface   : %d km2  " % S + '({:.0%})'.format(S/S0))
print("L moyenne : %.0f km   " % (S/((rn_h_strait-rn_h_ridge)/1E3) ) )
