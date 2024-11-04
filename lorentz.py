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

nx = int(paras[13])
ny = int(paras[14])
nz = int(paras[15])

nx= 128; ny = 128; nz = 128

nsnaps = int(paras[3])

remote_flag = int(paras[16])

xs = np.linspace(-12,12, nx+1)
ys = np.linspace(-12,12, ny+1)
zs = np.linspace(-24/nz,24, nz+1)

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

#data_sources = ['./Data_unstratified_64/', './Data_stratified_64/']
#data_sources = ['./Data/', './Data/']

#data_sources = ['/extra/tmp/trcn27/mf3d/001/','/extra/tmp/trcn27/mf3d/002/']
data_sources = ['/nobackup/trcn27/mf3d0/000/','/nobackup/trcn27/mf3d0/001/','/nobackup/trcn27/mf3d0/002/','/nobackup/trcn27/mf3d0/003/','/nobackup/trcn27/mf3d0/004/']


i = 0

if len(sys.argv) > 2:
    i = int(sys.argv[2])

snap_num = 0

cs = ['green', 'red', 'blue', 'orange', 'purple']

for plot_num in range(0, 120):
    ymax = 0.0; ymin = 1e6
    nulls = []
    for source_number in range(5):

        data_directory = data_sources[source_number]

        slice_index = ny//2
        i = plot_num
        wait = 0
        fname = '%s%04d.nc' % (data_directory, i)
        print('Examining plot', i, 'fname', fname)

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


        data.close()

        def current_test(bx,by,bz):
            jx = (bz[1:-1,1:,:] - bz[1:-1,:-1,:])/dy - (by[1:-1,:,1:] - by[1:-1,:,:-1])/dz
            jy =  (bx[:,1:-1,1:] - bx[:,1:-1,:-1])/dz - (bz[1:,1:-1,:] - bz[:-1,1:-1,:])/dx
            jz =  (by[1:,:,1:-1] - by[:-1,:,1:-1])/dx - (bx[:,1:,1:-1] - bx[:,:-1,1:-1])/dy

        def magfield(bx, by, bz):
            bx1 = 0.5*(bx[1:,slice_index,1:-1] + bx[:-1,slice_index,1:-1])
            by1 = 0.5*(by[1:-1,slice_index,1:-1] + by[1:-1,slice_index,1:-1])
            bz1 = 0.5*(bz[1:-1,slice_index,1:] + bz[1:-1,slice_index,:-1])
            return 0.5*(bx1**2 + by1**2+ bz1**2)


        current_test(bx, by, bz)

        def null_height(bz):
            #Determine(roughly) the height of the null point, up to the eruption
            zslice = bz[slice_index,slice_index,:]
            flip_index = np.where(np.sign(zslice[1:])*np.sign(zslice[:-1]) < 0)
            return 0.5*(zs[flip_index[-1]] + zs[flip_index[-1]+1])

        #Calculate Lorentz Force everywhere. Let's average to grid CENTRES because why not?
        jx = (bz[1:-1,1:,:] - bz[1:-1,:-1,:])/dy - (by[1:-1,:,1:] - by[1:-1,:,:-1])/dz
        jy =  (bx[:,1:-1,1:] - bx[:,1:-1,:-1])/dz - (bz[1:,1:-1,:] - bz[:-1,1:-1,:])/dx
        jz =  (by[1:,:,1:-1] - by[:-1,:,1:-1])/dx - (bx[:,1:,1:-1] - bx[:,:-1,1:-1])/dy

        #Average magnetic field to centres
        bx1 = 0.5*(bx[1:,1:-1,1:-1] + bx[:-1,1:-1,1:-1])
        by1 = 0.5*(by[1:-1,1:,1:-1] + by[1:-1,:-1,1:-1])
        bz1 = 0.5*(bz[1:-1,1:-1,1:] + bz[1:-1,1:-1,:-1])

        #Average current to centres
        jx1 = 0.25*(jx[:,1:,1:] + jx[:,1:,:-1] + jx[:,:-1,1:] + jx[:,:-1,:-1])
        jy1 = 0.25*(jy[1:,:,1:] + jy[1:,:,:-1] + jy[:-1,:,1:] + jy[:-1,:,:-1])
        jz1 = 0.25*(jz[1:,1:,:] + jz[1:,:-1,:] + jz[:-1,1:,:] + jz[:-1,:-1,:])


        lx1 = jy1*bz1 - jz1*by1
        ly1 = jz1*bx1 - jx1*bz1
        lz1 = jx1*by1 - jy1*bx1

        l2 = lx1**2 + ly1**2 + lz1**2

        #Run through vertical slices of these and output average values.
        #Can do Fourier transform if this makes sense, but will do so less than the flux rope sims

        lf_slice_avgs = np.zeros((nz))
        for k in range(nz):
            lf_slice_avgs[k] = np.sum(np.abs(l2[1:-1,1:-1,k]))/(np.size(l2[1:-1,1:-1,k]))

        plt.plot(zc[1:-1], lf_slice_avgs[1:-1], label = data_directory, c = cs[source_number])
        ymax = max(ymax, np.max(lf_slice_avgs[1:-1]))
        ymin = min(ymin, np.min(lf_slice_avgs[1:-1]))

        null = null_height(bz)
        nulls.append(null)
        #plt.plot([null, null], [ymin,ymax], linestyle = 'dashed', c = cs[source_number])

    for n in range(5):
        plt.plot([float(n), float(n)], [ymin,ymax], linestyle = 'dotted', c = cs[n])
        plt.plot([nulls[n], nulls[n]], [ymin,ymax], linestyle = 'dashed', c = cs[n])

    plt.plot([0.0, 0.0], [ymin,ymax], linestyle = 'solid', c = 'black')
    plt.xlabel('Height')
    plt.ylabel('Absolute Lorentz Force')
    plt.title('Output number %d' % plot_num)
    plt.legend()
    plt.savefig('./lorentz/%04d.png' % plot_num)
    plt.close()























