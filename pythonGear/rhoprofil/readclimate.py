# -*- coding: utf-8 -*-

import time, sys, os
import numpy as np
import gsw
import matplotlib.pyplot as plt

lat1=-2.25 ; lon1 = -29.
lat2=+2.4  ; lon2 = -10.

# quarter of a degree
# longitude ; latitude ; depth
idim=1440 ; jdim=720 ; kdim=33
N = idim*jdim*kdim

print("....TEMPERATURE")

fname= 'rhoprofil/t000AV_v2'
print("....reading %s" % fname)
f = open(fname, 'r+')
lines = [line for line in f.readlines()]
f.close()

# reorder data
tmptab=[]
for i in range(int(N/10)):
    if i%(N/100)==0:
        sys.stdout.write(u"\u001b[1000D" + "[%d/%d]" % (i,N))
        sys.stdout.flush()
    tmptab+=[float(lines[i][j*8:(j+1)*8]) for j in range(10)]

print("\n....reindexing data")
tarr = np.zeros((jdim,idim,kdim))
for k in range(kdim):
    for j in range(jdim):
        ip=(k+j)*idim ; ip1=(k+j+1)*idim
        tarr[j,:,k] = tmptab[ip:ip1]

print("....SALINITY")
fname= 'rhoprofil/s000AV_v2'
print("....reading %s" % fname)
f = open(fname, 'r+')
lines = [line for line in f.readlines()]
f.close()

# reorder data
tmptab=[]
for i in range(int(N/10)):
    if i%(N/100)==0:
        sys.stdout.write(u"\u001b[1000D" + "[%d/%d]" % (i,N))
        sys.stdout.flush()
    tmptab+=[float(lines[i][j*8:(j+1)*8]) for j in range(10)]

print("\n....reindexing data")
sarr = np.zeros((jdim,idim,kdim))
for k in range(kdim):
    for j in range(jdim):
        ip=(k+j)*idim ; ip1=(k+j+1)*idim
        sarr[j,:,k] = tmptab[ip:ip1]

print("......................................")

print("....building positions lattices")
nav_lon=np.zeros((idim)) ; nav_lat=np.zeros((jdim)) ; dx=1./4
for i in range(idim):
    nav_lon[i]= -90+i*dx + dx/2.
for j in range(jdim):
    nav_lat[j]= 0.+j*dx + dx/2.

# in situ temperature and salinity
ist = np.ma.masked_where(tarr[:] <= -99.9999,tarr[:])
iss = np.ma.masked_where(sarr[:] <= -99.9999,sarr[:])
print("....converting to density function")
rho = gsw.density.rho_t_exact(SA=iss, t=ist, p=4000.) - 1000.
