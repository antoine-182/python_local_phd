# -*- coding: utf-8 -*-

import netCDF4 as nc4
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
""" ********  FILES TO PLOT  *************
    *****************************************
    tenseur divrot - freeslip
    convergence 1°/4, /8, /16, /48
    0°, 45°
"""

dirt = "/Users/gm/Documents/nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXPREF/"
dirm = "/Users/gm/Documents/nemo/dev_r12527_Gurvan_ShallowWater/cfgs/AM98/EXPREF/"

ListA = [[4 , 0. , dirt+"een/AM98_1m_00010101_00041230_grid_T.nc",
                   dirm+"mesh_mask_0.nc"],
         [4 , 0. , dirt+"een_holl/AM98_1m_00010101_00041230_grid_T.nc",
                   dirm+"mesh_mask_0.nc"],
        [4 , 0. , dirt+"een_ren/AM98_1m_00010101_00041230_grid_T.nc",
                  dirm+"mesh_mask_0.nc"] ]


Listdt = [ListA]

""" *************************** PARAMETERS
"""
L=2000E3 # 2000 km
vmin = 200 ; vmax = 800

tickszer = [np.float64(i)*L/4. for i in range(5)]
cticks = np.arange(vmin, vmax+20., 5*20.)
# palette = plt.get_cmap('RdYlBu_r')
palette = plt.get_cmap('RdBu_r', 20)


save = 1 ; psave = "sshAMgrid.png" ; dpi = 200

NI = 3 ; NJ = 1 # nb of rows and cols

fig, ax = plt.subplots(NI,NJ, figsize=(NJ*3,NI*3), dpi = dpi, squeeze=False)

""" *************************** DATA
"""
for i in range(NI):
    for j in range(NJ):
        n,theta,pdtT,pmm = Listdt[j][i]
        dx = 100E3/np.float64(n)
        # reading
        print("... dataframe[%d,%d], n = %d" % (i,j,n))
        dt = nc4.Dataset(pdtT)
        mm = nc4.Dataset(pmm)
        tmask = mm.variables['tmask'][0,0]
        glamt = mm.variables['glamt'][0] ; gphit = mm.variables['gphit'][0]
        # grid
        lx = dx*np.cos((45+theta)*np.pi/180)/np.sqrt(2) ; ly = dx*np.sin((45+theta)*np.pi/180)/np.sqrt(2)
        nx,ny = np.shape(glamt)
        gridx = np.zeros((nx,ny)) ; gridy = np.zeros((nx,ny))
        gridx = glamt - lx ; gridy = gphit - ly
        # data
        ssh = dt.variables['sossheig'][-1,:,:] +500.
        data = np.ma.masked_where(tmask==0,ssh)
        # plot
        cf = ax[i][j].pcolormesh(gridx, gridy, data,
                            vmin=vmin, vmax=vmax, alpha = 0.9,
                            # levels = levels,
                            cmap = palette)
        # ax[i][j].title.set_text(titles[i][j])
        # contour
        c1= ax[i][j].contour(glamt,gphit, data,
                        levels = cticks,
                        vmin=vmin,vmax=vmax,
                        linewidths =0.3, colors=('k',),linestyles = "solid")
        ax[i][j].clabel(c1, fmt='%3.0f', colors='k', fontsize=8)


""" *************************** PRETTY
"""
for i in range(NI):
    for j in range(NJ):
        n,t,pdtT,pmm = Listdt[j][i]
        dx = 100E3/np.float64(n)
        limx = dx/2 ; limy = dx/2
        # print(tickszer)
        ax[i][j].set_xlim(-limx,L + limx)
        ax[i][j].set_xticks(tickszer)
        ax[i][j].set_ylim(-limy,L + limy)
        ax[i][j].set_yticks(tickszer)
        ax[i][j].set_xticklabels([]) ; ax[i][j].set_yticklabels([])
        #
        if (i==NI-1):
            ax[i][j].set_xticklabels(["%d" % (x/1E3) for x in tickszer])
            ax[i][j].set_xlabel("X (km)")
        #
        if (j==0):
            ax[i][j].set_yticklabels(["%d" % (x/1E3) for x in tickszer])
            ax[i][j].set_ylabel("Y (km)")
        #
        ax[i][j].patch.set_color('0.')
        ax[i][j].xaxis.set_minor_locator(AutoMinorLocator())
        ax[i][j].yaxis.set_minor_locator(AutoMinorLocator())
        ax[i][j].tick_params(axis = "y", which = 'both', width=1.5, labelsize = 10, pad = 5, left = True, right = True)
        ax[i][j].tick_params(axis = 'x', which = 'both', width=1.5, labelsize = 10, pad = 10, bottom = True, top = True)
        ax[i][j].tick_params(which='minor',length = 4)
        ax[i][j].tick_params(which='major',length = 6)
        ax[i][j].set_aspect(aspect='equal') # data coordinate 'equal'


fig.subplots_adjust(hspace = 0., wspace = 0.2, right=0.85)

""" *************************** COLORBAR
"""

for i in range(NI):
    pos= ax[i][-1].get_position()
    cax = fig.add_axes([pos.x1+0.03, pos.y0, 0.015, pos.y1-pos.y0])
    cbar = fig.colorbar(cf, cax=cax,
                        orientation = "vertical", extend = 'both')
    cbar.set_ticks(cticks)
    cbar.set_ticklabels(["%dm" % s for s in cticks])
    # cbar.set_label(r"Upper Layer Width - (m)")

""" *************************** SAVING
"""

if save :
    print("\nsaving : %s" % psave)
    fig.savefig(psave, dpi = dpi)
    plt.close()
    print("figure closed")
plt.show()
