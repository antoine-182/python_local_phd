# -*- coding: utf-8 -*-

import netCDF4 as nc4
import numpy as np
import time, sys

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
import types

dt = 10.
dx = 1E3 ; dz = 20. ; w_bvp = 1 # Transport W pénalisé 1(T) ou non (F)
print("dt = %d s ; dx = %d km ; dz = %d m ; w_bvp = %d" % (dt,dx,dz, w_bvp))

 ##########################################################################################

# pdout = "/Users/gm/Documents/nemo/release-4.0.1/tests/my_overflow/EXP_sanity/output.abort.nc"
# dout = nc4.Dataset(pdout)
# nav_lon = dout.variables['nav_lon'][1,:]
#
# vozocrtx = dout.variables['vozocrtx'][0,:,1,:] ; Cx = dt/dx
# vovecrtz = dout.variables['vovecrtz'][0,:,1,:] ; Cz = dt/dz
#
# rpot = dout.variables['rpot'][0,:,1,:]
# rpou = dout.variables['rpou'][0,:,1,:]
# rpow = dout.variables['rpow'][0,:,1,:]
#
# nK,nI = np.shape(vozocrtx)

##########################################################################################

# zWT = np.copy(vozocrtx)*0. ; zUT = np.copy(vozocrtx)*0.
# zWU = np.copy(vozocrtx)*0. ; zUU = np.copy(vozocrtx)*0.
#
# for jk in range(1,nK-1):
#     for ji in range(1,nI-1):
#         """ linéaire en T """
#         wp_kp =   0.5 * ( vovecrtz[jk+1,ji] + np.abs(vovecrtz[jk+1,ji]) )
#         wm_km = - 0.5 * ( vovecrtz[jk  ,ji] - np.abs(vovecrtz[jk,ji  ]) )
#         zWT[jk,ji] = Cz * np.max((rpow[jk+1,ji] * wp_kp / rpot[jk,ji], \
#                                   rpow[jk  ,ji] * wm_km / rpot[jk,ji]) )   #   CFL w -> T
#
#         up_im =   0.5 * ( vozocrtx[jk,ji-1] + np.abs(vozocrtx[jk,ji-1]) )
#         um_ip = - 0.5 * ( vozocrtx[jk,ji  ] - np.abs(vozocrtx[jk,ji  ]) )
#         zUT[jk,ji] = Cx * np.max((rpou[jk,ji-1] * up_im / rpot[jk,ji], \
#                                   rpou[jk,ji  ] * um_ip / rpot[jk,ji]) )   #   CFL u -> T
#         """ non-linéaire en U (rentrant dans la maille !)"""
#         wp = 0.5*(rpow[jk+1,ji]*vovecrtz[jk+1,ji]+rpow[jk+1,ji+1]*vovecrtz[jk+1,ji+1])
#         wp_kp =   0.5 * ( wp + np.abs(wp) )
#         #
#         wm = 0.5*(rpow[jk  ,ji]*vovecrtz[jk  ,ji]+rpow[jk  ,ji+1]*vovecrtz[jk  ,ji+1])
#         wm_km = - 0.5 * ( wm - np.abs(wm) )
#         #
#         zWU[jk,ji] = Cz * np.max((wp_kp / rpou[jk,ji], \
#                                   wm_km / rpou[jk,ji]) )                     #   CFL Wu -> U
#         #
#         up = 0.5*(rpou[jk,ji]*vozocrtx[jk,ji] + rpou[jk,ji-1]*vozocrtx[jk,ji-1])
#         up_im =   0.5 * ( up + np.abs(up) )
#         #
#         um = 0.5*(rpou[jk,ji+1]*vozocrtx[jk,ji+1] + rpou[jk,ji]*vozocrtx[jk,ji])
#         um_ip = - 0.5 * ( um - np.abs(um) )
#         #
#         zUU[jk,ji] = Cx * np.max((up_im / rpou[jk,ji], \
#                                   um_ip / rpou[jk,ji])  )                     #   CFL Uu -> U
#         # zUU[jk,ji] = Cx * up_im / rpou[jk,ji]                      #   CFL Uu -> U
#
# print("zWT.max() = %2.2f" % (zWT.max()))
# print("zUT.max() = %2.2f" % (zUT.max()))
#
# print("\nzWU.max() = %2.2f" % (zWU.max()))
# print("zUU.max() = %2.2f" % (zUU.max()))
#
# z2d = zWU
# print("\nz2d.max() = %2.2f" % (z2d.max()))
#
# fig, ax = plt.subplots()
# plt.imshow(z2d, cmap = plt.get_cmap('magma_r'))
# plt.xlim(10,70)
# # https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="2%", pad=0.05)
# cbar = plt.colorbar(cax=cax)
# cbar.set_label("cfl")
# plt.show()

##########################################################################################
# z2d = np.copy(rpou) * 0.
# tmask = np.copy(rpou) * 0. + 1.
# tmask[rpot<=1e-4]=0
#
# for jk in range(1,nK-1):
#     for ji in range(1,nI-1):
#         # c = 0.5 * ( rpou[jk,ji] + rpou[jk,ji-1] ) / rpou[jk,ji]  # --> non linéaire
#         #c = 0.5 * ( rpow[jk,ji] + rpow[jk,ji+1] ) / rpou[jk,ji]  # --> non linéaire
#         # c = rpow[jk,ji] / rpot[jk,ji]                            # --> linéaire
#         # c = rpou[jk,ji-1] / rpot[jk,ji]                            # --> linéaire
#         #
#         # c = 0.5 * ( rpot[jk,ji] + rpot[jk,ji-1] ) / rpot[jk,ji]  # --> non linéaire
#         #
#         # c =  tmask[jk,ji] * np.max( [rpou[jk,ji], rpou[jk,ji-1]] ) / rpot[jk,ji]  #
#         # c =  tmask[jk,ji] * rpou[jk,ji] / np.min( [ rpot[jk,ji], rpot[jk,ji+1] ] ) #
#         c =  np.max( [ 0.5*(rpou[jk,ji-1]+rpou[jk,ji]), 0.5*(rpou[jk,ji]+rpou[jk,ji+1]) ] ) / rpou[jk,ji]  #
#         # c =  np.max( [ 0.5*(rpow[jk,ji-1]+rpow[jk,ji]), 0.5*(rpow[jk,ji]+rpow[jk,ji+1]) ] ) / rpot[jk,ji]  #
#         z2d[jk,ji] = c
#
# # z2d = np.ma.array(z2d,mask=z2d<=1.)
# fig, ax = plt.subplots()
# plt.imshow(z2d, cmap = plt.get_cmap('magma_r'),vmin=1.)
# plt.xlim(10,70)
# # https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="2%", pad=0.05)
# cbar = plt.colorbar(cax=cax)
# cbar.set_label("max(phiu)/phit")
# plt.show()
##########################################################################################
dx = 1000. ; dz = 20. # km ; m
gridx = np.arange(0,200E3,dx) ; gridw = np.arange(0,2020,dz)
Ni = len(gridx) ; Nk = len(gridw)

def profilz(x) :
    return( 500. + 1500./2. * ( 1. + np.tanh( (x - 40.)/7. ) ) )

zht = np.copy(gridx)*0.
zht[:] = profilz(gridx/1E3)

def shapiro_2d(tab,n=1):
    z2d = np.copy(tab) ; (nz,nx) = np.shape(tab)
    for _ in range(n) :
        for k in range (nz) :
            for i in range (1,nx-1) :
                z2d[k,i] = 0.25 * tab[k,i-1] + 0.5 * tab[k,i] + 0.25 * tab[k,i+1]
        tab=z2d
        for k in range (1,nz-1) :
            for i in range (nx) :
                z2d[k,i] = 0.25 * tab[k-1,i] + 0.5 * tab[k,i] + 0.25 * tab[k+1,i]
        tab=z2d
    return(tab)


rn_abp = 1e-3

rpot = np.zeros((Nk,Ni))
for ji in range(Ni):
    rpot[gridw>zht[ji],ji ] = rn_abp
rpot[-1,:] = rn_abp

a = shapiro_2d(rpot)

x = np.arange(-10,10,0.5)
rpo0 = np.copy(x)
rpo0[a<=0] = 1.e-4 ; rpo0[a>0] = 1.
N = len(rpo0)

rpot=shapiro(rpo0,15) ; rpou=np.copy(rpot)
zp = np.copy(rpot)*0. ; zm = np.copy(rpot)*0. ; zw = np.copy(rpot)*0.
for i in range(1,N-1):
    rpou[i] = 0.5*(rpot[i]+rpot[i+1])

for i in range(1,N-1):
    zw[i] = 0.5*(rpot[i] + rpot[i+1])/rpou[i]

plt.plot(x,rpo0)
plt.plot(x,rpot)
plt.plot(x,zw,label=r"$\overline{\phi_w}/\phi_u$")
plt.show()

##########################################################################################

""" Tout dépend de la définition de Cu_adv
"""

# cc = np.arange(0.,10.,0.1)
#
# Cmax = 0.30 ; Cmin = 0.15
# Cut = 2*Cmax - Cmin
# Fcu = 4*Cmax*(Cmax-Cmin)
#
# fnemo = np.copy(cc)*0.
# for i in range(len(cc)):
#     if cc[i] <= Cmin:  # fully explicit
#         fnemo[i] = 0
#     elif cc[i] < Cut : # mixed implicit
#         z1d = (cc[i] - Cmin)**2
#         fnemo[i] =  z1d / ( Fcu + z1d )
#     else :
#         fnemo[i] = ( cc[i] - Cmax ) / cc[i]
#     fnemo[i] = np.min([1.,fnemo[i]])
#
# fsch = np.copy(cc)*0.  # 1-zcff
# for i in range(len(cc)):
#     if cc[i] <= Cmin:  # fully explicit
#         fsch[i] = 1.
#     elif cc[i] < Cut : # mixed implicit
#         z1d = (cc[i] - Cmin)**2
#         fsch[i] =  1. / ( 1. + z1d / Fcu )
#     else :
#         fsch[i] = Cmax / cc[i]
#     fsch[i] = 1 - np.min([1.,fsch[i]])
#
# plt.plot(cc,fnemo, label = "Nemo")
# plt.plot(cc,fsch , '.', label = "Shchepetkin")
# plt.vlines(Cmin,0,1)
# plt.vlines(Cut,0,1)
# plt.vlines(1,0,1)
# plt.xlim(0,cc[-1])
# plt.ylim(0,1)
# plt.legend()
# plt.xlabel("Courant Number")
# plt.ylabel("zcff")
# plt.grid()
# plt.show()
