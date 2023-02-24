# -*- coding: utf-8 -*-

import time, sys, os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy import interpolate

############

# Brasil Basin
ibb = 129 ; jbb = 68
# Serri Leone Abyssal Plane
isl = 138 ; jsl = 75

""" Build an idealised density profil across the romanche ridge
"""

pdt="ORCA2.init.nc"

dt = nc4.Dataset(pdt)
rho      = dt.variables['sigma4'  ][0,:,:,:]
vosaline = dt.variables['vosaline'][0,:,:,:]
votemper = dt.variables['votemper'][0,:,:,:]
gdepth   = dt.variables['gdepth'  ][0,:,:,:]

# nZ,nY,nX = np.shape(rho)
# ici le mask est là ou s=0 g/kg

# vertical reference profile of density (correction par rapport à rho_0 = 1028 kg/m3)
pp = 4000. * 1E-4 #  m to decibar
R00 = 4.6494977072e+01 ; R01 = -5.2099962525e+00 ; R02 =  2.2601900708e-01
R03 = 6.4326772569e-02 ; R04 =  1.5616995503e-02 ; R05 = -1.7243708991e-03
r0 = ( ( ( ( ( R05*pp + R04 )*pp + R03 )*pp + R02 )*pp + R01 )*pp + R00 )*pp
# r0 = 0.

zerBB   = vosaline[:,jbb,ibb]      ;  zerSL   = vosaline[:,jsl,isl]
rhoBB   = rho     [:,jbb,ibb] + r0 ;  rhoSL   = rho     [:,jsl,isl] + r0
depthBB = gdepth  [:,jbb,ibb]      ;  depthSL = gdepth  [:,jsl,isl]
toceBB  = votemper[:,jbb,ibb]      ;  toceSL  = votemper[:,jsl,isl]

rhoBB   = rhoBB  [zerBB != 0.]  ;  rhoSL   = rhoSL  [zerSL != 0.]
depthBB = depthBB[zerBB != 0.]  ;  depthSL = depthSL[zerSL != 0.]
toceBB  = toceBB [zerBB != 0.]  ;  toceSL  = toceSL [zerSL != 0.]

############
# h0 = 3180 #m
h0 = 3000. # rn_height

# h0 = 2200 # m
# depth[iBB] = 3250>h0 ce qui assure qu'il n'y aura pas de mélange
iBB    = (np.abs(depthBB - h0   )).argmin() ; iSL    = (np.abs(depthSL - h0   )).argmin()
iBB = 26 ; iSL = 26 # force for 3000m
iz1000 = (np.abs(depthBB - 1000.)).argmin() ; iz1000 = (np.abs(depthSL - 1000.)).argmin()

############
# before the 1000m

# from 1000m to h0 keep the same density -> define a function which goes between them
rho_up =           0.5 * (rhoBB  [iz1000-1:iz1000+2] + rhoSL  [iz1000-1:iz1000+2])
zhd    = np.log10( 0.5 * (depthBB[iz1000-1:iz1000+2] + depthSL[iz1000-1:iz1000+2]) )

rho_up  = np.append([         0.5 * (rhoBB  [0] + rhoSL  [0])], rho_up)
zhd     = np.append([np.log10(0.5 * (depthBB[0] + depthSL[0]) )], zhd   )   # ---> approche difficlement la densité au fond


fdeg, res, _, _,_ = np.polyfit(zhd, rho_up, deg = 3, full = True)
fup = np.poly1d(fdeg)

###########
# from 1000m to h0 keep the same density -> define a function which goes between them
rho_mid =           0.5 * (rhoBB  [iz1000:iBB+1] + rhoSL  [iz1000:iSL+1])
zhd     = np.log10( 0.5 * (depthBB[iz1000:iBB+1] + depthSL[iz1000:iSL+1]) )

# rho_mid =           0.5 * (rhoBB  [iz1000:] + rhoSL  [iz1000:])
# zhd     = np.log10( 0.5 * (depthBB[iz1000:] + depthSL[iz1000:]) )    # ---> donne sensiblement la même chose

# rho_mid = np.append([         0.5 * (rhoBB  [0] + rhoSL  [0])], rho_mid)
# zhd     = np.append([np.log10(0.5 * (depthBB[0] + depthSL[0]) )], zhd   )   # ---> approche difficlement la densité au fond

# rho_mid =           0.5 * (rhoBB  [:] + rhoSL  [:])
# zhd     = np.log10( 0.5 * (depthBB[:] + depthSL[:]) )    # ---> provoque un retournement

fdeg, res, _, _,_ = np.polyfit(zhd, rho_mid, deg = 3, full = True)
fmid  = np.poly1d(fdeg)
dfmid = np.poly1d.deriv(fmid)
############
dz = 10.
zdep = np.arange(h0,depthBB[iBB]+dz,dz)

h = depthBB[iBB] - h0
rho0 = fmid(np.log10(h0))  ;  dho0 = ( fmid(np.log10(zdep[1])) - fmid(np.log10(zdep[0])) ) / dz
rhoN = rhoBB[iBB]          ;  dhoN = (rhoBB[iBB+1]-rhoBB[iBB])/(depthBB[iBB+1]-depthBB[iBB])

fl = np.poly1d(np.array(( ( + (dhoN +   dho0)*h - 2*(rhoN-rho0) )/(h*h*h),
                          ( - (dhoN + 2*dho0)*h + 3*(rhoN-rho0) )/(h*h),
                               dho0,
                               rho0
                        )) )

h = depthSL[iSL] - h0
rho0 = fmid(np.log10(h0))  ;  dho0 = ( fmid(np.log10(zdep[1])) - fmid(np.log10(zdep[0])) ) / dz
rhoN = rhoSL[iSL]          ;  dhoN = (rhoSL[iSL+1]-rhoSL[iSL])/(depthSL[iSL+1]-depthSL[iSL])

fr = np.poly1d(np.array(( ( + (dhoN +   dho0)*h - 2*(rhoN-rho0) )/(h*h*h),
                          ( - (dhoN + 2*dho0)*h + 3*(rhoN-rho0) )/(h*h),
                               dho0,
                               rho0
                        )) )

############

dz = 10.
zhd_l = np.arange(depthBB[iBB],depthBB[-2]+dz,dz)

h = depthBB[-2] - depthBB[iBB]
rho0 = fl(depthBB[iBB]-h0)  ;  dho0 = ( fl(zhd_l[1]-h0) - fl(zhd_l[0]-h0) ) / dz
rhoN = rhoBB[-2]                   ;  dhoN = (rhoBB[-1]-rhoBB[-2])/(depthBB[-1]-depthBB[-2])

fdol = np.poly1d(np.array(( ( + (dhoN +   dho0)*h - 2*(rhoN-rho0) )/(h*h*h),
                            ( - (dhoN + 2*dho0)*h + 3*(rhoN-rho0) )/(h*h),
                                 dho0,
                                 rho0
                         )) )

zhd_r = np.arange(depthSL[iSL],depthSL[-2]+dz,dz)
h = depthSL[-2] - depthSL[iSL]
rho0 = fr(depthSL[iSL]-h0)  ;  dho0 = ( fr(zhd_r[1]-h0) - fr(zhd_r[0]-h0) ) / dz
rhoN = rhoSL[-2]       ;  dhoN = (rhoSL[-1]-rhoSL[-2])/(depthSL[-1]-depthSL[-2])

fdor = np.poly1d(np.array(( ( + (dhoN +   dho0)*h - 2*(rhoN-rho0) )/(h*h*h),
                            ( - (dhoN + 2*dho0)*h + 3*(rhoN-rho0) )/(h*h),
                                 dho0,
                                 rho0
                         )) )

##########
############
dz = 5.
refdepth = np.arange(depthBB[0], depthBB[-1]+2000., dz)

def rhd_left(z):
    if (z<=1000.):
        res = fup(np.log10(z))
    elif (z>1000. and z<=h0):
        res = fmid(np.log10(z))
    elif (z>h0 and z<depthBB[iBB]):
        res = fl(z-h0)
    elif (z>=depthBB[iBB]):
        res = fdol(z-zhd_l[0])
    return(res)

def rhd_right(z):
    if (z<=1000.):
        res = fup(np.log10(z))
    elif (z>1000. and z<=h0):
        res = fmid(np.log10(z))
    elif (z>h0 and z<depthBB[iBB]):
        res = fr(z-h0)
    elif (z>=depthBB[iBB]):
        res = fdor(z-zhd_r[0])
    return(res)


density_left = np.copy(refdepth)*0.
density_right = np.copy(refdepth)*0.
for jk in range(len(refdepth)):
    density_left[jk]=rhd_left(refdepth[jk])
    density_right[jk]=rhd_right(refdepth[jk])

#########
#########

# with open("denleft.txt", "w") as txt_file:
#     for jk in range(len(refdepth)):
#         line = ["%s" % refdepth[jk], "%s" % density_left[jk]]
#         txt_file.write(" ".join(line) + "\n") # works with any number of elements in a line
#
# with open("denright.txt", "w") as txt_file:
#     for jk in range(len(refdepth)):
#         line = ["%s" % refdepth[jk], "%s" % density_right[jk]]
#         txt_file.write(" ".join(line) + "\n") # works with any number of elements in a line


plt.plot(rhoBB      ,depthBB, 'o-', linewidth = 0., c = "royalblue")
plt.plot(rhoSL      ,depthSL, 'o-', linewidth = 0., c = "red")
#
plt.plot(fup (np.log10(depthSL))      ,depthSL, '--', c = "green")
plt.plot(fmid(np.log10(depthSL))      ,depthSL, '--', c = "black")
plt.plot(fl  (zdep-h0)         ,zdep, '--', c = "royalblue")
plt.plot(fr  (zdep-h0)         ,zdep, '--', c = "red")
plt.plot(fdol(zhd_l - zhd_l[0])      ,zhd_l, '--', c = "royalblue")
plt.plot(fdor(zhd_r - zhd_r[0])      ,zhd_r, '--', c = "red")

plt.plot(density_left      ,refdepth, '-', linewidth = 1., c = "royalblue")
plt.plot(density_right     ,refdepth, '-', linewidth = 1., c = "red")
#
plt.hlines(h0, 0, 50)
plt.gca().invert_yaxis()
plt.xlim((-0.1 + np.min([rhoBB[0],rhoSL[0]]), np.max([rhoBB[-1],rhoSL[-1]]) +0.1) )
plt.xlabel("potential density")
plt.ylabel("depth")
plt.grid()
plt.savefig("rho.png")
plt.show()
