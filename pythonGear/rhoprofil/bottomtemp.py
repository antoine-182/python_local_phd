import time, sys, os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4



pdt="ORCA2.init.nc"

dt = nc4.Dataset(pdt)
temp = dt.variables['votemper'][0,:,:,:]
nZ,nY,nX = np.shape(temp)
# ici le mask est là ou s=0 g/kg
vosaline = dt.variables['vosaline'][0,:,:,:]
data = np.ma.masked_where(vosaline == 0., temp)
############
