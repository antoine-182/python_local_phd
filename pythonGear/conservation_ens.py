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
import types
import numpy.ma as ma

""" *****************************************************************
"""

mean_vor = "mean4"

# ########################################################
#
dir  = "/Users/gm/Documents/nemo/dev_r12527_Gurvan_ShallowWater/cfgs/VORTEX/EXP_ens"

# name = "0_meanwet"
# dir_dt = dir+"/%s" % (name)
# pmm = dir_dt +"/mesh_mask.nc"
# mm = nc4.Dataset(pmm)
# tmask = mm.variables['tmask'][0][0,:,:]
# nY,nX = np.shape(tmask)
# # glamt = mm.variables['glamt'][0][:,:]
# fmask = np.copy(tmask)*0. ; ssfmask = np.copy(tmask)*0.
# r1_4 = np.float32(0.25) # ...
# for i in range(nX-1):
#     for j in range(nY-1):
#         #
#         fmask[j,i] = tmask[j+1,i] * tmask[j+1,i+1]  \
#                    * tmask[j  ,i] * tmask[j  ,i+1]
#         #
#         ssfmask[j,i] = np.max(  [  tmask[j+1,i] , tmask[j+1,i+1]  \
#                                  , tmask[j  ,i] , tmask[j  ,i+1]  ] )
# #
# pdt = dir_dt +"/VORTEX_1d_00010101_00010130_grid_T.nc"
# # pdu = dir_dt +"/VORTEX_1m_00010101_00041230_grid_U.nc"
# # pdv = dir_dt +"/VORTEX_1m_00010101_00041230_grid_V.nc"
# dtt = nc4.Dataset(pdt) # ; dtu = nc4.Dataset(pdu) ; dtv = nc4.Dataset(pdv)
# ht  = dtt.variables['ht'][:,:,:]
# nT,_,_ = np.shape(ht)
# hf = np.copy(ht)   *0.
# Af = np.copy(tmask)*0.
# #
# # if name.split('_')[-1] == 'mean4' :
# # for i in range(nX-1):
# #     for j in range(nY-1):
# #         # ht est déja maské
# #         hf[:,j,i] = r1_4 * (  ht[:,j+1,i] + ht[:,j+1,i+1] \
# #                             + ht[:,j  ,i] + ht[:,j  ,i+1])
# #         # 1 si un mouillé
# #         Af[j,i] = np.max(  [  tmask[j+1,i] , tmask[j+1,i+1]  \
# #                             , tmask[j  ,i] , tmask[j  ,i+1]  ] )
# # elif name.split('_')[-1] == 'meanwet' :
# for i in range(nX-1):
#     for j in range(nY-1):
#         sumtmask = np.sum(  [ tmask[j+1,i] , tmask[j+1,i+1]  \
#                             , tmask[j  ,i] , tmask[j  ,i+1]  ] )
#         #
#         hf[:,j,i] = (  ht[:,j+1,i] + ht[:,j+1,i+1] \
#                      + ht[:,j  ,i] + ht[:,j  ,i+1]) / sumtmask
#         #
#         Af[j,i] = sumtmask / 4
#
# tmp = np.sum(Af[:,:] * ssfmask[:,:] / ( hf[:,:,:] + 1. - ssfmask[:,:]) , axis=(1,2,))
#
# plt.plot(tmp)
# plt.show()


# time_sery[:nT,n] = tmp
# ma.masked_where(time_sery[:,n] == 0., time_sery[:,n])


testcase = ['0_mean4','0_meanwet','0_nl_mean4','0_nl_meanwet','45_mean4','45_meanwet']
# testcase = ['45_meanwet']
# testcase = ['mean4']
N = len(testcase)
time_sery = np.zeros((30,N))*np.nan

for n in range(N) :
    name = testcase[n]
    dir_dt = dir+"/%s" % (name)
    pmm = dir_dt +"/mesh_mask.nc"
    mm = nc4.Dataset(pmm)

    tmask = mm.variables['tmask'][0][0,:,:]
    nY,nX = np.shape(tmask)
    # glamt = mm.variables['glamt'][0][:,:]
    fmask = np.copy(tmask)*0. ; ssfmask = np.copy(tmask)*0.

    r1_4 = np.float32(0.25) # ...

    for i in range(nX-1):
        for j in range(nY-1):
            #
            fmask[j,i] = tmask[j+1,i] * tmask[j+1,i+1]  \
                       * tmask[j  ,i] * tmask[j  ,i+1]
            #
            ssfmask[j,i] = np.max(  [  tmask[j+1,i] , tmask[j+1,i+1]  \
                                     , tmask[j  ,i] , tmask[j  ,i+1]  ] )
    #

    pdt = dir_dt +"/VORTEX_1d_00010101_00010130_grid_T.nc"
    # pdu = dir_dt +"/VORTEX_1m_00010101_00041230_grid_U.nc"
    # pdv = dir_dt +"/VORTEX_1m_00010101_00041230_grid_V.nc"
    dtt = nc4.Dataset(pdt) #; dtu = nc4.Dataset(pdu) ; dtv = nc4.Dataset(pdv)
    ht  = dtt.variables['ht'][:,:,:]
    nT,_,_ = np.shape(ht)
    hf = np.copy(ht)   *0.
    Af = np.copy(tmask)*0.
    #
    # if name.split('_')[-1] == 'mean4' :
    # for i in range(nX-1):
    #     for j in range(nY-1):
    #         hf[:,j,i] = r1_4 * (  ht[:,j+1,i] + ht[:,j+1,i+1] \
    #                             + ht[:,j  ,i] + ht[:,j  ,i+1])
    #         #
    #         Af[j,i] = np.max(  [  tmask[j+1,i] , tmask[j+1,i+1]  \
    #                             , tmask[j  ,i] , tmask[j  ,i+1]  ] )
    # elif name.split('_')[-1] == 'meanwet' :
    for i in range(nX-1):
        for j in range(nY-1):
            sumtmask = np.sum(  [ tmask[j+1,i] , tmask[j+1,i+1]  \
                                , tmask[j  ,i] , tmask[j  ,i+1]  ] )
            #
            hf[:,j,i] = (  ht[:,j+1,i] + ht[:,j+1,i+1] \
                         + ht[:,j  ,i] + ht[:,j  ,i+1]) / sumtmask
            #
            Af[j,i] = sumtmask / 4.
    #
    tmp = np.sum(Af[:,:] * ssfmask[:,:] / ( hf[:,:,:] + 1. - ssfmask[:,:]) , axis=(1,2,))
    time_sery[:nT,n] = tmp
    # ma.masked_where(time_sery[:,n] == 0., time_sery[:,n])

fig, ax = plt.subplots()
for n in range(N):
    name = testcase[n]
    plt.plot(time_sery[:,n]/time_sery[0,n], label=name)
    # plt.plot(time_sery[:,n], label=name)
plt.legend()
ax.set_ylabel("vorticity tracer (passive)")
ax.set_xlabel("time (day)")

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", width=1.5, labelsize = 10, pad = 5, left = True)
ax.tick_params(axis = 'x', width=1.5, labelsize = 10, pad = 10, bottom = True)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)

plt.show()
