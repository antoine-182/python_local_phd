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

""" *****************************************************************
"""
# pgear convAM -n AM98_bvpK4_10y_1m_00010101_00101230_grid_T.nc AM98_bvpK8_10y_1m_00010101_00101230_grid_T.nc -m mesh_mask10y.nc -nspace 4 8

import argparse
parser = argparse.ArgumentParser(description = "To test the convergence of AM98")
parser.add_argument('-n','--netfile', nargs='+',
                    type=str, required=True,
                    help="Netcdf file path")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")

parser.add_argument("-nspace", type=float, nargs='+',
                    default=4,
                    help="grid spacing in km")
parser.add_argument("-d", type=int,
                    default=4,
                    help="lenght (number of cells -gc) of the penalisation")
parser.add_argument("-t","--theta", type=float,
                    default=0.,
                    help="rotation angle of the grid (default 0°)")
parser.add_argument("-s","--save", action="count", default=0,
                    help="save the current figure")
parser.add_argument("-ps","--psave", type=str,
                    default="zer",
                    help="name of the saved figure")
parser.add_argument("-nmean", type=int,
                    default=-17,
                    help="period for the final mean state if (>0)  \
                    the mean is computed from the last nmean frames, \
                    if (<0) the last nmean periods diagnosited from the peacks,\
                    if (=0) nothing is done.")
parser.add_argument("-m","--meshmask", type=str,
                    default="mesh_mask_025_0.nc",
                    help="meshmask associated to the file")
args = parser.parse_args()

if args.verbose:
    print("Running '{}'".format(__file__))
    print("netcdf used " + args.netfile)
nmean = args.nmean
listpdt = args.netfile
listnn = args.nspace
save=args.save
psave=args.psave+".png"
theta = float(args.theta)
pmm = args.meshmask
dx = 100E3 # 100 km
# theta = 0.

print(listpdt)
print(listnn)
print("nmean %d" % nmean)
""" ********  COARSEST RESOLUTION  *************
"""
# pmm  = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_bvp4/bvpK_s=1e-4/mesh_mask10y.nc"
nn = listnn[0]
mm = nc4.Dataset(pmm)
tmask = mm.variables['tmask'][0,0]
nI,nJ = np.shape(tmask)
glamt = mm.variables['glamt'][0]
gphit = mm.variables['gphit'][0]
Npoints = np.nansum(tmask)
Ndt = len(listpdt)
listLR = np.zeros((Ndt,nI,nJ))

""" ********  START REDUCTION  *************
"""
for ndt in range(Ndt):
    pdtT_hd=listpdt[ndt]
    dt_hd = nc4.Dataset(pdtT_hd)
    #
    nn_hd = listnn[ndt]
    #
    ssh_hd = dt_hd.variables['sossheig'][:,:,:]
    nT,ni,nj=np.shape(ssh_hd)
    glamt_hd = dt_hd.variables['nav_lon']
    gphit_hd = dt_hd.variables['nav_lat']

    print("\nDataframe %d/%d : n=%d" % (ndt+1,Ndt,nn_hd))

    """ ********  Mean Final State  *************
        *****************************************
    """
    # 256 (period de 15.1 mois)
    if nmean>0:
        # mean on a multiple of the period
        ssh_hd = np.mean(ssh_hd[(nT-nmean):,:,:], axis=0)
    elif nmean<0:
        # variance
        Ep = np.zeros((nT))
        for t in range(nT):
            Ep[t]=np.sum(ssh_hd[t])
        # peacks
        # dataraw = np.flip(Ep)[:400]  # not compatible with numpy 1.7
        dataraw = Ep[::-1][:400]
        data = dataraw - np.mean(dataraw)
        index = np.diff(np.sign(np.diff(data)))
        extrema = np.where(index>0.)[0] +1 # minima
        # take the last nmean period
        t0=nT-extrema[0] ; t = nT-extrema[nmean]
        ssh_hd=np.mean(ssh_hd[t:t0,:,:], axis=0)
    else:
        ssh_hd = ssh_hd[-1,:,:]
    plt.imsave('pre%d.pdf' % (ndt),ssh_hd,vmin=-200, vmax=200, dpi=400)
    #
    """ ********  HD -> Coarse Resolution  *************
        ************************************************
    """
    if nn!=nn_hd:
        LR = np.copy(tmask)*0.
        alpha = 0.5/nn ; beta = (nn/nn_hd)**2
        jj = 0 ; JJ = 0 ; counter = 0
        while (JJ<nJ and jj<nj) :
            ii = 0 ; II = 0
            if   (glamt_hd[ii,jj] < glamt[II,JJ] - dx*alpha) :
                jj+=1   # derrière JJ
            elif (glamt_hd[ii,jj] > glamt[II,JJ] + dx*alpha) :
                JJ+=1   # derrière jj
            while (II<nI and ii<ni):
                if   (gphit_hd[ii,jj] < gphit[II,JJ] - dx*alpha) :
                    ii+=1   # derrière II
                elif (gphit_hd[ii,jj] > gphit[II,JJ] + dx*alpha) :
                    II+=1   # derrière ii
                else :
                    counter+=1
                    sys.stdout.write(u"\u001b[1000D" + "processing [%8d/%8d]" % (counter,ni*nj))
                    sys.stdout.flush()
                    LR[II,JJ]+= beta * ssh_hd[ii,jj]
                    ii+=1
            jj+=1
        # end while
    else:
        LR = ssh_hd
    plt.imsave('%d.pdf' % (ndt),LR,vmin=-200, vmax=200, dpi=400)
    # SAVE
    listLR[ndt,:,:] = LR

""" ********  CALCULS DES ERREURS  *************
"""
nl2 = np.zeros((Ndt-1)) ; nl1 = np.zeros((Ndt-1)) ; ninf = np.zeros((Ndt-1))

for ndt in range(Ndt-1):
    data = listLR[-1]-listLR[ndt]
    # data = listLR[ndt]
    nl2 [ndt]=np.sqrt(np.sum(data*data))/Npoints
    nl1 [ndt]=np.sum(np.abs(data))/Npoints
    ninf[ndt]=np.nanmax(np.abs(data))
xnl2 = [t for t in listnn[:-1]]

with open("nerror.txt", "w") as txt_file:
    for ndt in range(Ndt-1):
        line = ["%s" % xnl2[ndt], "%s" % nl2[ndt], "%s" % nl1[ndt], "%s" % ninf[ndt]]
        txt_file.write(" ".join(line) + "\n") # works with any number of elements in a line

fig, ax = plt.subplots(dpi=200)
ax.plot(xnl2, nl2,
        marker="o", color="royalblue",
        alpha=1., linewidth=.7, markersize=5.)
ax.grid(which='major', linestyle='-', linewidth='0.3', color='black')
ax.grid(which='minor', linestyle=':', linewidth='0.3', color='black')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis = "y", which = 'both', width=1., labelsize = 10, pad = 5)
ax.tick_params(axis = 'x', which = 'both', width=1., labelsize = 10, pad = 10)
ax.tick_params(which='minor',length = 4)
ax.tick_params(which='major',length = 6)
ax.set_ylabel("L2 error")
ax.set_xlabel("Resolution")
plt.xscale('log')
plt.yscale('log')
fig.subplots_adjust(top = 0.8, bottom=0.15, hspace = 0.02)

if save :
    print("saving : %s" % psave)
    fig.savefig(psave)
else :
    plt.show()


#
# """ *****************************************************************
# """
#
#
# """ ********  COARSEST RESOLUTION  *************
# """
#
# pdtT = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_bvp4/bvpK_s=1e-4/AM98_bvpK4_10y_1m_00010101_00101230_grid_T.nc"
# pmm  = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_bvp4/bvpK_s=1e-4/mesh_mask10y.nc"
# nn   = 4
#
# dt = nc4.Dataset(pdtT)
# mm = nc4.Dataset(pmm)
#
# ssh = dt.variables['sossheig'][:,:,:]
# tmask = mm.variables['tmask'][0,0]
# nT,nI,nJ = np.shape(ssh)
# glamt = mm.variables['glamt'][0]
# gphit = mm.variables['gphit'][0]
#
# """ ********  HIGHEST RESOLUTION  *************
#     *****************************************
# """
#
# pdtT_hd = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_bvp4/bvpK8_s=1e-4/AM98_bvpK8_10y_1m_00010101_00101230_grid_T.nc"
# # pmm_hd  = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_bvp4/bvpK8_s=1e-4/mesh_mask.nc"
# nn_hd   = 8
#
# # pdtT_hd = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_bvp4/bvpK_s=1e-4/AM98_bvpK4_10y_1m_00010101_00101230_grid_T.nc"
# # nn_hd   = 4
#
# dt_hd = nc4.Dataset(pdtT_hd)
#
# ssh_hd = dt_hd.variables['sossheig'][:,:,:]
# _,ni,nj = np.shape(ssh_hd)
# glamt_hd = dt_hd.variables['nav_lon']
# gphit_hd = dt_hd.variables['nav_lat']
#
# """ ********  Mean Final State  *************
#     *****************************************
# """
# # cf clean.py la moyenne sur un multiple entier de 14.22 mois fonctionne (71)
# ssh    = np.mean(ssh   [(nT-71):,:,:], axis=0)
# ssh_hd = np.mean(ssh_hd[(nT-71):,:,:], axis=0)
#
# """ ********  HD -> Coarse Resolution  *************
#     ************************************************
# """
# LR = ssh*0.
# alpha = 0.5/nn ; beta = (nn/nn_hd)**2
# jj = 0 ; JJ = 0 ; counter = 0
# while (JJ<nJ and jj<nj) :
#     ii = 0 ; II = 0
#     if   (glamt_hd[ii,jj] < glamt[II,JJ] - dx*alpha) :
#         jj+=1   # derrière JJ
#     elif (glamt_hd[ii,jj] > glamt[II,JJ] + dx*alpha) :
#         JJ+=1   # derrière jj
#     while (II<nI and ii<ni):
#         if   (gphit_hd[ii,jj] < gphit[II,JJ] - dx*alpha) :
#             ii+=1   # derrière II
#         elif (gphit_hd[ii,jj] > gphit[II,JJ] + dx*alpha) :
#             II+=1   # derrière ii
#         else :
#             counter+=1
#             sys.stdout.write(u"\u001b[1000D" + "processing [%8d/%8d]" % (counter,ni*nj))
#             sys.stdout.flush()
#             LR[II,JJ]+= beta * ssh_hd[ii,jj]
#             ii+=1
#     jj+=1
#
# # print("counter : %d/%d" % (counter,ni*nj))
#
# """ ********  PLOT FIELD + DIFF  ************
#     *****************************************
# """
# vmin = -100 ; vmax = 100 ; alph = 0.9
# dx/=nn
# # grid centré sur les corners bas gauche
# lx = dx*np.cos((45+theta)*np.pi/180)/np.sqrt(2)
# ly = dx*np.sin((45+theta)*np.pi/180)/np.sqrt(2)
# # ballec cotes droit/haut - ils seront cropés
# nx = np.shape(glamt)[0] ; ny = np.shape(glamt)[1]
#
# gridx = np.zeros((nx,ny))
# gridx = glamt - lx
# gridy = np.zeros((nx,ny))
# gridy = gphit - ly
#
# data = LR-ssh
# data = np.ma.masked_where(tmask==0,data)
# titlezer = 'zer\n'
# psave = "zer.png"
#
# palette = plt.get_cmap('RdYlBu_r')
# # palette = plt.get_cmap('seismic')
# """ figure bordelum """
# fig, ax = plt.subplots(dpi=200)
# im = ax.pcolormesh(gridx, gridy, data,
#                     vmin=vmin, vmax=vmax, alpha = alph,
#                     # levels = levels,
#                     cmap = palette)
# cbar = plt.colorbar(im)
# cbar.set_label(r"Upper Layer Width - (m)")
#
# # adapté pour tout les angles
# limx = dx/2 ; limy = dx/2 ; d = 4
# bvpgapx = 2*lx * float(d) ; bvpgapy = 2*ly * float(d)
#
#
# ax.set_xlim(-limx - bvpgapx ,2000000 + limx + bvpgapx )
# ax.set_xticks([0,500000,1000000,1500000,2000000])
# ax.set_xticklabels(["0","500","1000","1500","2000"])
# ax.set_xlabel("X (km)")
#
# ax.set_ylim(-limy - bvpgapy ,2000000 + limy + bvpgapy )
# ax.set_yticks([0,500000,1000000,1500000,2000000])
# ax.set_yticklabels(["0","500","1000","1500","2000"])
# ax.set_ylabel("Y (km)")
#
# if d > 0 :
#     print("add true coastlines")
#     ax.axvline(x= 0     , ymin = 0, ymax = 2000E3, linestyle = '-', color = 'red', lw = 1.)
#     ax.axvline(x= 2000E3, ymin = 0, ymax = 2000E3, linestyle = '-', color = 'red', lw = 1.)
#     ax.axhline(y= 0     , xmin = 0, xmax = 2000E3, linestyle = '-', color = 'red', lw = 1.)
#     ax.axhline(y= 2000E3, xmin = 0, xmax = 2000E3, linestyle = '-', color = 'red', lw = 1.)
#
# titlezer += "min = %3.1f m   max = %3.1f m   " % (np.min(data), np.max(data)) \
#          +  "$\Delta$ = %3.1f m" % (np.max(data) - np.min(data))
# ax.set_title(titlezer,
#               fontsize = 10, y = 1.02)
# ##################################################
#
# ax.patch.set_color('0.')
# ax.xaxis.set_minor_locator(AutoMinorLocator())
# ax.yaxis.set_minor_locator(AutoMinorLocator())
# ax.tick_params(axis = "y", which = 'both', width=1.5, labelsize = 10, pad = 5, left = True, right = True)
# ax.tick_params(axis = 'x', which = 'both', width=1.5, labelsize = 10, pad = 10, bottom = True, top = True)
# ax.tick_params(which='minor',length = 4)
# ax.tick_params(which='major',length = 6)
# ax.set_aspect(aspect=1) # data coordinate 'equal'
#
# plt.show()
#
# """ ******** L2 error  *************
#     ********************************
# """
#
# nl2=np.sum(data*data)
