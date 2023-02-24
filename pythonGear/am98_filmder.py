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

import argparse
parser = argparse.ArgumentParser(description = "Plot the vorticiy field for Deremble's (2016) testcase")
parser.add_argument("netfile", type=str,
                    help="Netcdf file path")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")

parser.add_argument("-n", type=int,
                    default=4,
                    help="grid spacing in km")
parser.add_argument("-o", type=int,
                    default=0,
                    help="grid spacing in km of obstacles")
parser.add_argument("-nhd", type=int,
                    default=0,
                    help="grid spacing in km of hd solution (for reducing)")
parser.add_argument("-d", type=int,
                    default=0,
                    help="lenght (number of cells -gc) of the penalisation")
parser.add_argument("-t","--theta", type=float,
                    default=0.,
                    help="rotation angle of the grid (default 0°)")
parser.add_argument("-L","--domsize", type=float,
                    default=2000.,
                    help="daomaine size (default 2000km)")
parser.add_argument("-s","--save", action="count", default=0,
                    help="save the current figure")
parser.add_argument("-f","--film", action="count", default=0,
                    help="make a film")
parser.add_argument("-sta","--stack", action="count", default=0,
                    help="adapt the grid to the porosity field")

parser.add_argument("-min","--minimum", type=float, default = 200,
                    help="minimal value of the colorbar")
parser.add_argument("-max","--maximum", type=float, default = 800,
                    help="maximal value of the colorbar")
parser.add_argument("-mm","--minmax", action="count", default=0,
                    help="add the moment the gap tends toward its final value")
parser.add_argument("-a","--alpha", type=float,
                    default=.9,
                    help="transparency figure")
parser.add_argument("-duree","--duration", type=float,
                    default=20.,
                    help="lenght of the movie (in s)")

parser.add_argument("-m","--meshmask", type=str,
                    default="mesh_mask_025_0.nc",
                    help="meshmask associated to the file")

parser.add_argument("-nmean", type=int,
                    default=0,
                    help="period for the final mean state if (>0)  \
                    the mean is computed from the last nmean frames, \
                    if (<0) the last nmean periods diagnosited from the peacks,\
                    if (=0) nothing is done.")

args = parser.parse_args()

if args.verbose:
    print("Running '{}'".format(__file__))
    print("netcdf used " + args.netfile)

""" *****************************************************************
"""
theta = float(args.theta)
dx = 100E3 # 100km
dx/= float(args.n)
if args.o == 0:
    dxo = 100E3/float(args.n)
else :
    dxo = 100E3/float(args.o)
dy=dx
L = args.domsize*1E3
duree = args.duration
# dx =  float(args.dx) ; dy =  float(args.dx) # 25 km
# dx*=1000;dy*=1000
pdtT = args.netfile
mxmn = args.minmax
alph = args.alpha
print(args.meshmask)
# Result in, if there is no min/max declared, =None and pcolormesh handle this
vmin = args.minimum ; vmax = args.maximum
nmean = args.nmean
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

# pdtT = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_ref4/freeslip0/AM98_ref_4_freeslip_0_40y_5d_00010101_00401230_grid_T.nc"
# pmm  = "../nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXP_ref4/freeslip0/mesh_mask.nc"
# dx = 25 ; dy =  25 # 25 km
# dx*=1000;dy*=1000
# theta = 0.

print("nmean %d" % nmean)

""" *****************************************************************
"""

# ssh
dt = nc4.Dataset(pdtT)

try:
    ssh = dt.variables['sossheig'][:,:,:]
except:
    ssh = dt.variables['ht'][:,:,:]
nT,_,_ = np.shape(ssh)
# nT=100

mm = nc4.Dataset(pmm)
tmask = mm.variables['tmask'][0,0]
# glamt, gphit
glamt = mm.variables['glamt'][0]
gphit = mm.variables['gphit'][0]

# grid centré sur les corners bas gauche
lx = dx*np.cos((45+theta)*np.pi/180)/np.sqrt(2)
ly = dx*np.sin((45+theta)*np.pi/180)/np.sqrt(2)
# ballec cotes droit/haut - ils seront cropés
nx = np.shape(glamt)[0] ; ny = np.shape(glamt)[1]

""" Grille régulière/stackée """
if args.stack==0:
    gridx = np.zeros((nx,ny))
    gridx = glamt - lx
    gridy = np.zeros((nx,ny))
    gridy = gphit - ly
elif args.stack==1:
    gridx = np.zeros((nx,ny))
    gridx = glamt - lx
    gridy = np.zeros((nx,ny))
    gridy = gphit - ly

# adapté pour tout les angles
limx = dxo/2 ; limy = dxo/2
# limx = 0. ; limy = 0.
# limx = dx ; limy = dx
# limx = dx*np.cos((45+theta)*np.pi/180)
# limy = dx*np.sin((45+theta)*np.pi/180)
bvpgapx = 2*lx * float(args.d) ; bvpgapy = 2*ly * float(args.d)

""" ********  Mean Final State  *************
    *****************************************
"""
# 256 (period de 15.1 mois)
if nmean>0:
    # mean on a multiple of the period
    data = np.mean(ssh[(nT-nmean):,:,:], axis=0)
    #
    titlezer  = '%s - meaned over last %d frames \n'%(pdtT,nmean)
    psave = "ssh.png"
elif nmean<0:
    # variance
    Ep = np.zeros((nT))
    for t in range(nT):
        Ep[t]=np.sum(ssh[t])
    # peacks
    # dataraw = np.flip(Ep)[:400]  # not compatible with numpy 1.7
    dataraw = Ep[::-1][:400]
    data = dataraw - np.mean(dataraw)
    index = np.diff(np.sign(np.diff(data)))
    extrema = np.where(index>0.)[0] +1 # minima
    # take the last nmean period
    t0=nT-extrema[0] ; t = nT-extrema[-nmean]
    data=np.mean(ssh[t:t0,:,:], axis=0)
    #
    titlezer  = '%s - last %d periods \n'%(pdtT,-nmean)
    psave = "ssh%d.png" % (-nmean)
else:
    if args.film :
        data = ssh[0,:,:]
        titlezer  = '%s - %02d/%02d \n'%(pdtT,0+1,nT)
        psave = "ssh.mp4"
    else :
        data = ssh[-1,:,:]
        titlezer = '%s\n'% pdtT
        psave = "ssh.png"


""" ********  HD -> Coarse Resolution  *************
    ************************************************
"""
nn=float(args.n) ; nn_hd = float(args.nhd)
if nn_hd>nn:
    glamt_hd = dt.variables['nav_lon']
    gphit_hd = dt.variables['nav_lat']
    LR = np.copy(tmask)*0.
    nI,nJ = np.shape(tmask) ; _,ni,nj=np.shape(ssh)
    #
    alpha = 0.5 ; beta = (nn/nn_hd)**2
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
                LR[II,JJ]+= beta * data[ii,jj] # matrice de masse pour etre sur (nécessaire à la rotation)
                ii+=1
        jj+=1
    # end while
    data = LR

########################################################
data=data+500
data = np.ma.masked_where(tmask==0,data)

# if mxmn :
# t = np.arange(nT)
# minssh = np.array([np.nanmin(ssh[tt,:,:]) for tt in t] ) ; maxssh = np.array([np.nanmax(ssh[tt,:,:]) for tt in t] )
# delta = maxssh-minssh
# d1  = np.where(np.abs(delta-delta[-1]) <1)[0] ; d01 = np.where(np.abs(delta-delta[-1]) <0.1)[0]
# ad1 = np.argwhere(np.diff(d1)!=1)[-1][0]+1   ; ad01 = np.argwhere(np.diff(d01)!=1)[-1][0]+1
# titlezer += "(%d - %d / %d)\n" % (d1[ad1],d01[ad01],nT)

########################################################

cticks = np.arange(vmin, vmax+20., 5*20.)
# palette = plt.get_cmap('RdYlBu_r')
palette = plt.get_cmap('RdBu_r', 30)
# palette = plt.get_cmap('seismic')

""" figure bordelum """
fig, ax = plt.subplots(dpi=200)

im = ax.pcolormesh(gridx, gridy, data,
                    vmin=vmin, vmax=vmax, alpha = alph,
                    # levels = levels,
                    cmap = palette)
cbar = plt.colorbar(im)
if args.film==0 :
    c1 = ax.contour(glamt,gphit, data,
                    levels = cticks,
                    vmin=vmin,vmax=vmax,
                    linewidths =0.3, colors=('k',),linestyles = "solid")
    ax.clabel(c1, fmt='%3.0f', colors='k', fontsize=10)

cbar.set_label(r"Upper Layer Width - (m)")
cbar.set_ticks(cticks)
cbar.set_ticklabels(["%dm" % s for s in cticks])

# tickszer = [0,500000,1000000,1500000,2000000]
tickszer = [np.float64(i)*L/4. for i in range(5)]
# print(tickszer)
ax.set_xlim(-limx - bvpgapx ,L + limx + bvpgapx )
ax.set_xticks(tickszer)
ax.set_xticklabels(["%d" % (x/1E3) for x in tickszer])
ax.set_xlabel("X (km)")

ax.set_ylim(-limy - bvpgapy ,L + limy + bvpgapy )
ax.set_yticks(tickszer)
ax.set_yticklabels(["%d" % (x/1E3) for x in tickszer])
ax.set_ylabel("Y (km)")

if args.d > 0 :
    print("zer")
    ax.axvline(x= 0     , ymin = 0, ymax = L, linestyle = '-', color = 'red', lw = 1.)
    ax.axvline(x= L, ymin = 0, ymax = L, linestyle = '-', color = 'red', lw = 1.)
    ax.axhline(y= 0     , xmin = 0, xmax = L, linestyle = '-', color = 'red', lw = 1.)
    ax.axhline(y= L, xmin = 0, xmax = L, linestyle = '-', color = 'red', lw = 1.)

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
        data = 500.+ssh[n,:,:]
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
    writer = animation.writers['ffmpeg'](fps=nT/duree)
    ani.save(psave,writer=writer,dpi=200)
    print("\nsaving : %s" % psave)
else :
    if args.save :
        print("\nsaving : %s" % psave)
        fig.savefig(psave, dpi = 200)
        plt.close()
    plt.show()
# plt.show()
