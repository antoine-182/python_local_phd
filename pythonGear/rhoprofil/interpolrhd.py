# -*- coding: utf-8 -*-

import time, sys, os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy import interpolate

############

# Brasil Basin
# i1 = 129 ; j1 = 68 ; N = 7
# Serri Leone Abyssal Plane
i1 = 138 ; j1 = 75 ; N = 9

""" Search the polynomial fit of a given density profil
"""

pdt="ORCA2.init.nc"


dt = nc4.Dataset(pdt)
rho = dt.variables['sigma4'][0,:,j1,i1]
# nZ,nY,nX = np.shape(rho)
# ici le mask est là ou s=0 g/kg
vosaline = dt.variables['vosaline'][0,:,j1,i1]
# rho = np.ma.masked_where(vosaline == 0., rho)
gdepth = dt.variables['gdepth'][0,:,j1,i1]
############

rhd   = rho   [vosaline != 0.]
gdept = gdepth[vosaline != 0.]

zhd = np.log10(gdept)
gdeptnew = np.arange(gdept[0], gdept[-1], 20)
palette = plt.get_cmap('jet',N)

# sans préconditionnement les min oscillent, fonction
# f = interpolate.interp1d(gdept, rhd, kind = 'cubic')
for iin in range(N):
    n = iin+2
    if n%2==0:
        continue
    # fdeg, res, _, _,_ = np.polyfit(gdept, rhd, deg = n, full = True)
    # f = np.poly1d(fdeg)

    fdeg, res, _, _,_ = np.polyfit(zhd, rhd, deg = n, full = True)
    f = np.poly1d(fdeg)

    rho_inter = f(np.log10(gdeptnew))   # use interpolation function returned by `interp1d`
    plt.plot(rho_inter,gdeptnew, '-', label="n=%d" % n, color=palette(iin))
plt.plot(rhd      ,gdept   , 'o',linewidth = 0., c = "black")
# #########
# # plt.plot(rhd      ,gdept   , 'o-', c = "black")
# # plt.plot(rho_inter,gdeptnew, '-', c = "red" )
plt.gca().invert_yaxis()
plt.legend()
plt.show()

n = len(fdeg)
for i in range(n):
    # hisghest power first
    print("a_%2d = " % (n-i) + "{0:+E}".format(fdeg[i]) )


# plt.plot(BB      ,gdept, 'o-', c = "red")
# plt.plot(SL      ,gdept, 'o-', c = "royalblue")
# plt.gca().invert_yaxis()
# plt.xlabel("sigma4")
# plt.ylabel("depth")

""" Fonctions de log10(z) ou z est la profondeur en m
"""
# Brasil Basin
# a_ 0 = +1.019328E+00
# a_ 1 = -1.721174E+01
# a_ 2 = +1.204360E+02
# a_ 3 = -4.497664E+02
# a_ 4 = +9.599694E+02
# a_ 5 = -1.158413E+03
# a_ 6 = +7.240736E+02
# a_ 7 = -1.550587E+02

# SLAP
# a_ 0 = +6.618686E-01
# a_ 1 = -1.536698E+01
# a_ 2 = +1.548979E+02
# a_ 3 = -8.869686E+02
# a_ 4 = +3.167744E+03
# a_ 5 = -7.284099E+03
# a_ 6 = +1.072372E+04
# a_ 7 = -9.680916E+03
# a_ 8 = +4.827086E+03
# a_ 9 = -9.820019E+02
