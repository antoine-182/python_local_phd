
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

""" ********** Détermine le profil ****************
"""

dir = "/Users/gm/Documents/nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_bvp4/convref/"

# bord = "West" ; section = 1800E3
bord = "North" ; section = 900E3
L = 200E3 # km
nmean = 0
Listdt = [[0. ,  dir+"AM98_ref_4_fs_0_40y_1m_00010101_00401230_grid_U.nc",  \
                 dir+"AM98_ref_4_fs_0_40y_1m_00010101_00401230_grid_V.nc"], \
          [0. ,  dir+"AM98_ref_16_fs_0_40y_1m_00010101_00401230_grid_U.nc",  \
                 dir+"AM98_ref_16_fs_0_40y_1m_00010101_00401230_grid_V.nc"], \
          [45. , dir+"AM98_ref_4_fs_45_40y_1m_00010101_00401230_grid_U.nc", \
                 dir+"AM98_ref_4_fs_45_40y_1m_00010101_00401230_grid_V.nc"],  \
          [45. , dir+"AM98_ref_16_fs_45_40y_1m_00010101_00401230_grid_U.nc", \
                 dir+"AM98_ref_16_fs_45_40y_1m_00010101_00401230_grid_V.nc"],
          [45. , dir+"AM98_ref_domain_1m_00010101_00051230_grid_U.nc", \
                 dir+"AM98_ref_domain_1m_00010101_00051230_grid_V.nc"] ]
colorlist=["royalblue","royalblue","tab:red","tab:red","black"]
stylelist = ["o-","s--","o-","s--","v-"]

# profil Bord Ouest
# bord = "West" ; section = 1650E3
# bord = "North" ; section = 200E3
# L = 500E3 # km
# nmean = 3
# Listdt = [[0. ,  dir+"AM98_ref_4_ns_0_40y_1m_00010101_00401230_grid_U.nc",  \
#                  dir+"AM98_ref_4_ns_0_40y_1m_00010101_00401230_grid_V.nc"], \
#           [0. ,  dir+"AM98_ref_16_ns_0_40y_1m_00010101_00401230_grid_U.nc",  \
#                  dir+"AM98_ref_16_ns_0_40y_1m_00010101_00401230_grid_V.nc"], \
#           [45. , dir+"AM98_ref_4_noslip_45_5d_00010101_00051230_grid_U.nc", \
#                  dir+"AM98_ref_4_noslip_45_5d_00010101_00051230_grid_V.nc"],  \
#           [45. , dir+"AM98_ref_16_noslip_45_5d_00010101_00051230_grid_U.nc", \
#                  dir+"AM98_ref_16_noslip_45_5d_00010101_00051230_grid_V.nc"]]
# colorlist=["royalblue","royalblue","tab:red","tab:red"]
# stylelist = ["o-","s--","o-","s--"]

# bord = "West" ; section = 1800E3
# bord = "North" ; section = 500E3
# L = 200E3 # km
# nmean = 0
# Listdt = [[0. ,  dir+"AM98_ref_4_freeslip_nu2_0_40y_1m_00010101_00401230_grid_U.nc",  \
#                  dir+"AM98_ref_4_freeslip_nu2_0_40y_1m_00010101_00401230_grid_V.nc"], \
#           [0. ,  dir+"AM98_ref_8_freeslip_nu2_0_40y_1m_00010101_00401230_grid_U.nc",  \
#                  dir+"AM98_ref_8_freeslip_nu2_0_40y_1m_00010101_00401230_grid_V.nc"], \
#           [45. , dir+"AM98_ref_4_freeslip_nu2_45_40y_1m_00010101_00401230_grid_U.nc", \
#                  dir+"AM98_ref_4_freeslip_nu2_45_40y_1m_00010101_00401230_grid_V.nc"],  \
#           [45. , dir+"AM98_ref_8_freeslip_nu2_45_40y_1m_00010101_00401230_grid_U.nc", \
#                  dir+"AM98_ref_8_freeslip_nu2_45_40y_1m_00010101_00401230_grid_V.nc"]]
# colorlist=["royalblue","royalblue","tab:red","tab:red"]
# stylelist = ["o-","s--","o-","s--"]


psave = "%s%d.png" % (bord,section/1E3) ; save = 1


# Ni = 1
Ni=len(Listdt)

XY = np.zeros((Ni,300)) ; PU = np.zeros((Ni,300)) ; PT = np.zeros((Ni,300)) ; meta = np.zeros((Ni,4))
""" ***********************************************
"""
for dti in range(Ni):
    dataframe = Listdt[dti]
    # lecture
    angle = dataframe[0]
    cosangle = 1/np.cos(angle*np.pi/180) ; sinangle = 1/np.sin(angle*np.pi/180)
    # cosangle=1. ; sinangle = 1.
    dtU = nc4.Dataset(dataframe[1])
    dtV = nc4.Dataset(dataframe[2])

    uut = dtU.variables['ssu'][:,:,:]  ;  vvt = dtV.variables['ssv'][:,:,:]
    hu = dtU.variables['hu' ][:,:,:]  ;  hv = dtV.variables['hv' ][:,:,:]

    glamu = dtU.variables['nav_lon'] ; glamv = dtV.variables['nav_lon'] # longitude (x)
    gphiu = dtU.variables['nav_lat'] ; gphiv = dtV.variables['nav_lat'] # latitude  (y)
    ni,nj = np.shape(glamu)

    # nmean
    if nmean>0:
        nT,_,_ = np.shape(hu)
        Ep = np.zeros((nT))
        for t in range(nT):
            Ep[t]=np.sum(hu[t])
        # peacks
        # dataraw = np.flip(Ep)[:400]  # not compatible with numpy 1.7
        dataraw = Ep[::-1][:400]
        data = dataraw - np.mean(dataraw)
        index = np.diff(np.sign(np.diff(data)))
        extrema = np.where(index>0.)[0] +1 # minima
        # take the last nmean period
        t0=nT-extrema[0] ; t = nT-extrema[nmean]
        uu=np.mean(uut[t:t0,:,:], axis=0) ; vv=np.mean(vvt[t:t0,:,:], axis=0)
    else :
        uu = uut[-1] ; vv = vvt[-1]
    # #
    dl = np.diff(glamu[1,:])[0]      # résolution normale à la côte côte
    dx = dl/np.cos(angle*np.pi/180)  # e1

    print("angle : %02d°    dl = %2.1fkm     dx = %2.1fkm" % (angle,dl/1E3,dx/1E3))
    # Calcul
    if bord == "West":   # J ?
        if angle <0.:
            nn=np.int64(np.round(L/dx))
            # Determiner (ii,jj) XY
            j0 = np.where(np.abs(glamu[1,:]-0.     )<dl/2.)[0][0]+1
            i0 = np.where(np.abs(gphiv[:,1]-section)<dl/2.)[0][0]
            # Extraire le profil
            xy = glamv[i0,(j0-1):(j0+nn)]   #   ok pour angle = 0
            pu = vv   [i0,(j0-1):(j0+nn)]
            pt = vv   [i0,(j0-1):(j0+nn)]*hv[i0,(j0-1):(j0+nn)]
        else:
            # cherche un v mortel
            ii,jj = np.where( (np.abs(gphiv[:]-section)<dl/2.) \
                             *(       gphiv[:]<2000E3+dl)      \
                             *(       gphiv[:]>0.-dl))                 # si /dl j'ai la bande complète ca doit etre dl/2 resolution
            tmp1=np.zeros(len(ii)) ; tmpu1=np.zeros(len(ii)) ; tmpt=np.zeros(len(ii))
            for _ in range(len(ii)):
                tmp1[_]      =glamv[ii[_],jj[_]]
                tmpu1[_]     =   vv[ii[_],jj[_]]# /dl rattrape les u
            # cherche un u mortel
            ii,jj = np.where( (np.abs(gphiu[:]-section)<dl/2.) \
                             *(       gphiu[:]<2000E3+dl)      \
                             *(       gphiu[:]>0.-dl))                 # si /dl j'ai la bande complète ca doit etre dl/2 resolution
            tmp2=np.zeros(len(ii)) ; tmpu2=np.zeros(len(ii)) ; tmpt=np.zeros(len(ii))
            for _ in range(len(ii)):
                tmp2[_]      =glamu[ii[_],jj[_]]
                tmpu2[_]     =   uu[ii[_],jj[_]]# /dl rattrape les u
            # concatène
            tmp=np.zeros((len(tmp1)+len(tmp2))) ; tmpu = np.zeros((len(tmp1)+len(tmp2)))
            tmp [:(len(tmp1))] = tmp1  ; tmp [(len(tmp1)):]=tmp2
            tmpu[:(len(tmp1))] = tmpu1 ; tmpu[(len(tmp1)):]=tmpu2
            #sort
            index = np.argsort(tmp)
            tmp=tmp[index] ; tmpu=tmpu[index]
            # (2025,2000-L/1E3 )
            i0=np.where( (tmp<L)*(tmp>-dx*2) )[0]
            xy = tmp [i0]   ;   pu = tmpu[i0]   #;   pt = tmpt[i0]
            #
        meta[dti]=[angle,dl,dx,len(xy)]
        print(meta[dti])
        print("xy",xy)
        print("pu",pu)
        # print("pt",pt)
    elif bord=="North":
        if angle <0.:
            nn=np.int64(np.round(L/dx))
            # Determiner (ii,jj) XY
            j0 = np.where(np.abs(glamu[1,:]-section )<dl)[0][0]
            i0 = np.where(np.abs(gphiv[:,1]-2000E3 )<dl)[0][0]
            # Extraire le profil
            xy = gphiu[(i0-nn):(i0+2),j0]   #   ok pour angle = 0
            pu = uu   [(i0-nn):(i0+2),j0]
            # pt = vv   [i0,(j0-1):(j0+nn)]*hv[i0,(j0-1):(j0+nn)]
        else:
            # cherche un v mortel
            ii,jj = np.where( (np.abs(glamv[:]-section)<dl/2.) \
                             *(       glamv[:]<2000E3+dl)      \
                             *(       glamv[:]>0.-dl))                 # si /dl j'ai la bande complète ca doit etre dl/2 resolution
            tmp1=np.zeros(len(ii)) ; tmpu1=np.zeros(len(ii)) ; tmpt=np.zeros(len(ii))
            for _ in range(len(ii)):
                tmp1[_]      =gphiv[ii[_],jj[_]]
                tmpu1[_]     =   -vv[ii[_],jj[_]]# /dl rattrape les u
            # cherche un u mortel
            ii,jj = np.where( (np.abs(glamu[:]-section)<dl/2.) \
                             *(       glamu[:]<2000E3+dl)      \
                             *(       glamu[:]>0.-dl))
            tmp2=np.zeros(len(ii)) ; tmpu2=np.zeros(len(ii)) ; tmpt=np.zeros(len(ii))
            for _ in range(len(ii)):
                tmp2[_]      =gphiu[ii[_],jj[_]]
                tmpu2[_]     =   uu[ii[_],jj[_]]# /dl rattrape les u
            # concatène
            tmp=np.zeros((len(tmp1)+len(tmp2))) ; tmpu = np.zeros((len(tmp1)+len(tmp2)))
            tmp [:(len(tmp1))] = tmp1  ; tmp [(len(tmp1)):]=tmp2
            tmpu[:(len(tmp1))] = tmpu1 ; tmpu[(len(tmp1)):]=tmpu2
            #sort
            index = np.argsort(tmp)
            tmp=tmp[index] ; tmpu=tmpu[index]
            #
            i0=np.where( (tmp<2000E3+dx*2)*(tmp>2000E3-L) )[0]
            xy = tmp [i0]   ;   pu = tmpu[i0]   #;   pt = tmpt[i0]
            #
        meta[dti]=[angle,dl,dx,len(xy)]
        print(meta[dti])
        print("xy",xy)
        print("pu",pu)
    else:
        print("Error (bord) : boudary not declared")

    # remplissage array
    for _ in range(len(xy)):
        XY[dti,_] = xy[_]/1E3   ;   PU[dti,_] = pu[_]  # ; PT[dti,_] = pt[_]

# PU[0,:]*=0.7 # ATTENTION, entre 0° et 45° il y une projection
# PU[2,:]*=0.7
""" ***********************************************
    A 45°, le premier zéros est la facette maské contre la côte
"""
print("plot...")
fig, ax = plt.subplots(dpi=200)
for dti in range(Ni):
    dataframe = Listdt[dti]
    #
    angle,dl,dx,nn = meta[dti]
    print("angle : %02d°    dl = %2.1fkm     dx = %2.1fkm" % (angle,dl/1E3,dx/1E3))
    #
    mask=np.zeros(300)
    mask[:int(nn)]+=1   #  1 in profil
    #
    data=np.ma.masked_where(mask<1,PU[dti])
    ax.plot(XY[dti], data,stylelist[dti], color = colorlist[dti],
            label ="%02d°, dx=%02.2fkm, dl=%02.2fkm" % (angle,dx/1E3,dl/1E3),
            alpha=1., linewidth=.7, markersize=3.)
    data=np.ma.masked_where(mask<1,PT[dti])
    # ax.plot(XY[dti], data, "o--", color = color[dti],
    #         label ="_no_legend_",
    #         alpha=1., linewidth=.7, markersize=3.)
ax.grid(which='major', linestyle='-', linewidth='0.3', color='black')
ax.grid(which='minor', linestyle=':', linewidth='0.3', color='black')
# ax.set_ylim(-0.2,1.5 )
if bord=="West":
    ax.set_xlim(-25,L/1E3 )
    plt.vlines(0.,ymin=-0.2,ymax=2.,color="black")
elif bord=="North":
    ax.set_xlim(2025,2000-L/1E3 )
    plt.vlines(2000,ymin=-0.2,ymax=2.,color="black")

ax.xaxis.set_major_locator(MultipleLocator(25.  ))
ax.xaxis.set_minor_locator(MultipleLocator(25./4))
ax.yaxis.set_minor_locator(AutoMinorLocator())

ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
ax.set_ylabel("Tangential speed (m/s)")
ax.set_xlabel("Distance ashore (km)")
ax.set_title("Speed Profile against %s coast at y=%4dkm" % (bord,section/1E3), fontsize = 10, y = 1.02)
plt.legend()

""" Insert section in plot
"""
left, bottom, width, height = [.65, 0.35, 0.2, 0.25]
ax_new = fig.add_axes([left, bottom, width, height])
data = dtU.variables['hu'][-1,:,:]
data = np.ma.masked_where(data<1.,data)
ax_new.pcolormesh(glamu,gphiu,data,
              vmin=200, vmax=800, alpha = 0.9,
              cmap = plt.get_cmap('RdYlBu_r') )
if bord=="West":
    plt.hlines(section,xmin=xy[0],xmax=xy[-1],color="k")
elif bord=="North":
    plt.vlines(section,ymin=xy[0],ymax=xy[-1],color="k")
# plt.hlines(xyc[0],xmin=glamu[ii,jj],xmax=glamu[ii,jj],color="red")
ax_new.set_xlim(-dl ,2000000 + dl )
ax_new.set_xticks([0,500000,1000000,1500000,2000000])
ax_new.set_xticklabels(["0","500","1000","1500","2000"])
# ax_new.set_xlabel("X (km)")
ax_new.set_ylim(-dl ,2000000 + dl )
ax_new.set_yticks([0,500000,1000000,1500000,2000000])
ax_new.set_yticklabels(["0","500","1000","1500","2000"])
# ax_new.set_ylabel("Y (km)")

ax_new.patch.set_color('0.')
ax_new.xaxis.set_minor_locator(AutoMinorLocator())
ax_new.yaxis.set_minor_locator(AutoMinorLocator())
ax_new.tick_params(axis = "y", which = 'both', width=.5, labelsize = 0, pad = 5, left = True, right = True)
ax_new.tick_params(axis = 'x', which = 'both', width=.5, labelsize = 0, pad = 10, bottom = True, top = True)
ax_new.tick_params(which='minor',length = 2)
ax_new.tick_params(which='major',length = 3)
ax_new.set_aspect(aspect=1) # data coordinate 'equal'

if save :
    print("saving : %s" % psave)
    fig.savefig(psave)
    plt.close()
else :
    plt.show()
