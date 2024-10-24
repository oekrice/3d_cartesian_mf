#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: trcn27
"""

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import matplotlib
from scipy.io import netcdf_file
from fltrace import trace_fieldlines

#matplotlib.rcParams['text.usetex'] = True


if len(sys.argv) > 1:
    run = int(sys.argv[1])
else:
    run = 0

paras = np.loadtxt('parameters/variables%03d.txt' % run)

fig_width = 15#1.0*(513.11743/72)

nx = int(paras[1])
ny = int(paras[2])
nz = int(paras[3])

nsnaps = int(paras[5])

remote_flag = int(paras[22])

xs = np.linspace(-12,12, nx+1)
ys = np.linspace(-12,12, ny+1)
zs = np.linspace(0,24, nz+1)

xc = 0.5*(xs[1:] + xs[:-1])
yc = 0.5*(ys[1:] + ys[:-1])
zc = 0.5*(zs[1:] + zs[:-1])

dx = xs[1] - xs[0]
dy = ys[1] - ys[0]
dz = zs[1] - zs[0]

xc = 0.5*(xs[1:] + xs[:-1])
yc = 0.5*(ys[1:] + ys[:-1])
zc = 0.5*(zs[1:] + zs[:-1])

dx = xs[1] - xs[0]
dy = ys[1] - ys[0]
dz = zs[1] - zs[0]


class Grid():
    def __init__(self):
        self.x0 = xs[0]; self.x1 = xs[-1]
        self.y0 = ys[0]; self.y1 = ys[-1]
        self.z0 = zs[0]; self.z1 = zs[-1]
        self.nx = nx ; self.ny = ny; self.nz = nz

i = 0

if len(sys.argv) > 2:
    i = int(sys.argv[2])

for run in range(4):
#for plot_num in range(0,nsnaps,1):
    plot_num = 9
    if remote_flag:
        data_directory = '/nobackup/trcn27/mf3d0/%03d/' % run
    else:
        data_directory = '/extra/tmp/trcn27/mf3d/%03d/' % run

    slice_index = ny//2
    i = plot_num
    wait = 0
    fname = '%s%04d.nc' % (data_directory, i)
    print('Looking at file', i, 'fname', fname)

    try:
        data = netcdf_file(fname, 'r', mmap=False)
        print('File', fname, 'found')

    except:
        print('File', fname, 'not found')
        continue

    bx = np.zeros((nx+1,ny+2,nz+2))
    by = np.zeros((nx+2,ny+1,nz+2))
    bz = np.zeros((nx+2,ny+2,nz+1))

    bx[:,1:-1,1:-1] = np.swapaxes(data.variables['bx'][:],0,2)
    by[1:-1,:,1:-1] = np.swapaxes(data.variables['by'][:],0,2)
    bz[1:-1,1:-1,:] = np.swapaxes(data.variables['bz'][:],0,2)


    jx = np.zeros((nx+2,ny+1,nz+1))
    jy = np.zeros((nx+1,ny+2,nz+1))
    jz = np.zeros((nx+1,ny+1,nz+2))

    jx[1:-1,:,:] = np.swapaxes(data.variables['jx'][:],0,2)
    jy[:,1:-1,:] = np.swapaxes(data.variables['jy'][:],0,2)
    jz[:,:,1:-1] = np.swapaxes(data.variables['jz'][:],0,2)

    ex = np.zeros((nx+2,ny+1,nz+1))
    ey = np.zeros((nx+1,ny+2,nz+1))
    ez = np.zeros((nx+1,ny+1,nz+2))

    ex[1:-1,:,:] = np.swapaxes(data.variables['ex'][:],0,2)
    ey[:,1:-1,:] = np.swapaxes(data.variables['ey'][:],0,2)
    ez[:,:,1:-1] = np.swapaxes(data.variables['ez'][:],0,2)

    data.close()

    if run == 0:
        ex_reference = ex[1:-1,:,:]
        ey_reference = ex[:,1:-1,:]
        ez_reference = ex[:,:,1:-1]

    else:
        print('Test number', run)
        print('ex', np.allclose(ex_reference, ex[1:-1,:,:]))
        print('ey', np.allclose(ex_reference, ex[1:-1,:,:]))
        print('ez', np.allclose(ex_reference, ex[1:-1,:,:]))




